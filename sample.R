####################################################################################################
################ This script demonstrates an example for using our algorithm #######################
####################################################################################################
rm(list = ls(all = TRUE))
###################################################
library("glmnet")
library("matrixcalc")
library("fields")
library("vars")
library("MTS")
library("mvtnorm")
library("xtable")
library("lattice")
library("ggplot2")
library("doParallel")
library("factorcpt")

source("fista_LS.R")
source("functions_SBDetection.R")

###################################################
######## Generating dataset and parameters ########
###################################################
T <- 300 # number of time points
k <- 20  # number of time series
r <- 2
brk <- c(floor(T/3), floor(2*T/3), T+1)
m <- length(brk)
p.t <- 1; ## the true AR order
ratio <- 4;

### Generating low-rank component (for simulation we use fixed low-rank component loaded from file)
set.seed(10000)
L.true <- matrix(0, k, k)
basis.L <- matrix(rnorm(k*r), nrow = k, ncol = r)
coef.mat <- matrix(rnorm(k*r), nrow = r, ncol = k)
L.true <- basis.L %*% coef.mat
L.info <- norm(L.true, "m")

S.info <- L.info * ratio
sparse.true = phi.true <- matrix(0, nrow = k, ncol = m*k)
S.true <- vector("list", m)
for(i in 1:m){
    S.true[[i]] <- matrix(0, k, k)
    for(j in 1:(k-1)){
        S.true[[i]][j,j+1] <- 1
    }
}
S.true[[1]] <- S.true[[1]] * (-S.info)
S.true[[2]] <- S.true[[2]] * (S.info)
S.true[[3]] <- S.true[[3]] * (-S.info)

tmp.phi <- vector("list", m)
tmp.eigen <- c()
for(j in 1:m){
    tmp.phi[[j]] <- L.true + S.true[[j]]
    tmp.eigen <- c(max(abs(eigen(tmp.phi[[j]])$values)), tmp.eigen)
}
cnst <- max(tmp.eigen)
L.true <- L.true * 0.9 / cnst
for(j in 1:m){
    S.true[[j]] <- S.true[[j]] * 0.9 / cnst
    phi.true[,((j-1)*k+1):(j*k)] <- L.true + S.true[[j]]
    sparse.true[,((j-1)*k+1):(j*k)] <- S.true[[j]]
}

phi.full <- phi.true;
print(plot.matrix(abs(phi.full), m))

#########################################################
#################### Algorithm start ####################
#########################################################
### General paramters 
tol.sel <- 0.15
tol <- 2*10^(-2)                ### tolerances 
step.size <- 2*10^(-4)          ### step size 
max.iteration <- 100            ### max number of iteration for the LASSO solution
p <- p.t                        ### the selected AR order
sig <- matrix(0, k, k)
max.rep <- 1                    ### the number of total replications
e.sigma <- as.matrix(0.01*diag(k))

mu.sel.vec <- seq(2.5, 5, 5)
lambda.sel <- 0.5

lambda.1.cv <- seq(0.1, 0.5, 0.05)  ### The first tuning parameter
cv.length <- 10
cv.index <- seq(1+cv.length, T-1, cv.length) ### cv index

### Create result list to restore first step result
first.result = L.new.lst <- vector("list", max.rep)
for(i in 1:max.rep){
    ### Generating dataset
    set.seed(123456*i)
    try = var.sim.break(T, arlags=seq(1, p.t, 1), malags = NULL, phi = phi.full, sigma = e.sigma, brk = brk)
    data <- try$series
    data <- as.matrix(data)
    
    ################################################################################
    ############## First Step: Select candidate change points via L+S ##############
    ################################################################################
    L.old <- 0*diag(k)
    X.sel <- data[1:(T-1),]
    Y.sel <- data[2:T,]
    
    bic.seq <- rep(0, length(mu.sel.vec))
    L.new.list <- list()
    #### we employ 1-D grid search to find tuning parmaeter for low-rank component
    for(ll in 1:length(mu.sel.vec)){
        mu.sel <- mu.sel.vec[ll]
        fit.lr <- fista.nuclear(X.sel, Y.sel, lambda = mu.sel, k, niter = 100, backtracking = TRUE, t(L.true))
        L.new.list[[ll]] <- t(fit.lr$phi.hat)
        bic.seq[ll] <- f.func(L.new.list[[ll]], X.sel, Y.sel) + nuclear.pen(L.new.list[[ll]], mu.sel)
    }
    idx <- which.min(bic.seq)
    L.new <- L.new.list[[idx]]
    mu.sel <- mu.sel.vec[idx]
    qr(L.new)$rank
    
    if(norm(L.new - L.old, "F") < tol.sel){
        ### removing estimated low-rank effects from original data
        dat.sel <- matrix(0, T, k)
        dat.sel[1,] <- data[1,]
        X.tmp <- data[1:(T-1),]
        Y.tmp <- data[2:T,]
        dat.sel[2:T,] <- Y.tmp - X.tmp %*% t(L.new)
        
        ### Apply sparse fused-lasso to estimate change point candidates by using the removed low-rank data
        temp.first <- first.step.cv.new("LASSO", dat.sel, weight = NULL, lambda.1.cv, p, max.iteration = max.iteration, 
                                        tol = tol, step.size = (1/2)*10^(-4), cv.index)
        fisrt.brk.points <- temp.first$brk.points
        print(fisrt.brk.points)
        pts <- fisrt.brk.points
    }
    
    rep <- 0
    while(norm(L.new - L.old, "F") > tol.sel){
        rep <- rep + 1
        L.old <- L.new
        
        ### removing estimated low-rank effects from original data
        dat.sel <- matrix(0, T, k)
        dat.sel[1,] <- data[1,]
        X.tmp <- data[1:(T-1),]
        Y.tmp <- data[2:T,]
        dat.sel[2:T,] <- Y.tmp - X.tmp %*% t(L.new)
        
        ### Apply sparse fused-lasso to estimate change point candidates by using the removed low-rank data
        temp.first <- first.step.cv.new("LASSO", dat.sel, weight = NULL, lambda.1.cv, p, max.iteration = max.iteration, 
                                        tol = tol, step.size = (1/2)*10^(-4), cv.index)
        fisrt.brk.points <- temp.first$brk.points
        print(fisrt.brk.points)
        pts <- fisrt.brk.points
        
        ### separate data according to estimated time change points, and estimate sparse components
        m.hat <- length(pts)
        pts.est <- rep(0, m.hat+2)
        pts.est[1] <- 0
        pts.est[m.hat+2] <- T
        pts.est[2:(m.hat+1)] <- pts
        
        S.hat <- vector("list", m.hat+1)
        data.seg <- vector("list", m.hat+1)
        data.tmp <- vector("list", m.hat+1)
        for(jj in 1:(m.hat+1)){
            data.seg[[jj]] <- as.matrix(dat.sel[((pts.est[jj]+1):pts.est[jj+1]),]) ### removed low-rank effect data
            data.tmp[[jj]] <- as.matrix(data[((pts.est[jj]+1):pts.est[jj+1]),])    ### original data copy for estimated segments
        }
        
        ### separately estimate sparse components by using removed low-rank effect data
        X.sp = Y.sp <- vector("list", m.hat+1)
        for(jj in 1:(m.hat+1)){
            len <- dim(data.seg[[jj]])[1]
            X.sp[[jj]] <- data.tmp[[jj]][1:(len-1),]
            Y.sp[[jj]] <- data.seg[[jj]][2:len,]
            fit.S <- fista.sparse(X.sp[[jj]], Y.sp[[jj]], lambda = lambda.sel, k, niter = 40, backtracking = TRUE, 0*diag(k))
            S.hat[[jj]] <- t(fit.S$phi.hat)
        }
        
        ### Now remove sparse components from original data in each segments
        data.lr = X.lr = Y.lr <- vector("list", m.hat+1)
        for(jj in 1:(m.hat+1)){
            len <- dim(data.seg[[jj]])[1]
            data.lr[[jj]] <- matrix(0, dim(data.tmp[[jj]])[1], dim(data.tmp[[jj]])[2])
            data.lr[[jj]][1,] <- data.tmp[[jj]][1,]
            X.lr[[jj]] <- data.tmp[[jj]][1:(len-1),]
            Y.lr[[jj]] <- data.tmp[[jj]][2:len,]
            data.lr[[jj]][-1,] <- Y.lr[[jj]] - X.lr[[jj]] %*% t(S.hat[[jj]])
        }
        
        ### combine all removed sparse component data, and re-estimate low-rank matrix
        data.sp <- c()
        for(jj in 1:(m.hat+1)){
            data.sp <- rbind(data.sp, as.matrix(data.lr[[jj]]))
        }
        
        X.sel <- data[1:(T-1),]        ### design matrix should be from original data
        Y.sel <- data.sp[2:T,]         ### response matrix should be from removed effect data
        bic.seq <- rep(0, length(mu.sel.vec))
        L.new.list <- list()
        for(ll in 1:length(mu.sel.vec)){
            mu.sel <- mu.sel.vec[ll]
            fit.lr <- fista.nuclear(X.sel, Y.sel, lambda = mu.sel, k, niter = 100, backtracking = TRUE, t(L.true))
            L.new.list[[ll]] <- t(fit.lr$phi.hat)
            bic.seq[ll] <- f.func(L.new.list[[ll]], X.sel, Y.sel) + nuclear.pen(L.new.list[[ll]], mu.sel)
        }
        idx <- which.min(bic.seq)
        L.new <- L.new.list[[idx]]
        print(c(rep, norm(L.new - L.old, "F"), qr(L.new)$rank))
    }
    print(paste("Iters:", i))
    first.result[[i]] <- pts
    L.new.lst[[i]] <- L.new
}

############################################################################
############## Second Step: Screening candidates to get final ##############
############################################################################
second.result <- vector("list", max.rep)
for(i in 1:max.rep){
    ### Generating dataset
    set.seed(123456*i)
    try = var.sim.break(T, arlags=seq(1, p.t, 1), malags = NULL, phi = phi.full, sigma = e.sigma, brk = brk)
    data <- try$series
    data <- as.matrix(data)
    
    L.new <- L.new.lst[[i]]
    dat.screen <- matrix(0, T, k); dat.screen[1,] <- data[1,]
    X.tmp <- data[1:(T-1),]
    Y.tmp <- data[2:T,]
    dat.screen[-1,] <- Y.tmp - X.tmp %*% t(L.new)
    n <- T - 1
    
    ### the penalty term in the information criterion: the higher omega, the smaller number of break points selected.
    omega <- (1/3.5)*(((log(n))^1)*log(k))
    ### the second tuning parameter. Default number seems to be working for many simulation and real data examples!
    lambda.2 <- (1/1)*(log(n)*log(k))/n     
    
    ### employ the second step screening step for selecting final change points.
    temp <- second.step(dat.screen, lambda = lambda.2, p, max.iteration = 200, tol = tol, step.size = 10^(-3), first.result[[i]], omega)
    final.brk.points <- temp$pts  ### final break points selected
    
    ### plotting the data with selected break points
    MTSplot(data)
    abline(v = final.brk.points)
    pts.final <- final.brk.points
    second.result[[i]] <- pts.final
}
