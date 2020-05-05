rm(list = ls(all = TRUE))
############### Pre-processing #####################
library(imager)
library(plyr)

dir <- "~/surveillance/JPEGS"   
filenames <- list.files(dir, pattern = ".jpg", full.names = TRUE)

T <- length(filenames)
k <- 32*24
data <- matrix(0, nrow = T, ncol = k)

for(i in 1:length(filenames)){
    img <- load.image(filenames[i])
    gray.img <- grayscale(img)
    thmb <- resize(gray.img, round(width(img)/12), round(height(img)/12))
    dat <- thmb[,,1,]
    vdat <- as.vector(dat)
    data[i,] <- vdat
}

####################################################################################################
library("igraph")
library("ggplot2")
library("glmnet")
library("matrixcalc")
library("fields")
library("vars")
library("MTS")
library("mvtnorm")
library("xtable")
library("lattice")
library("doParallel")
source("fista_LS.R")
source("functions_SBDetection.R")

T <- dim(data)[1]
k <- dim(data)[2]
X.sel <- as.matrix(data[1:(T-1),])
Y.sel <- as.matrix(data[2:T,])

### Parameters settings 
tol.sel <- 0.5
tol <- 2*10^(-2)                ### tolerance 
step.size <- 2*10^(-4)          ### step size 
p = p.t <- 1                    ### the selected AR order

### Tuning parameter for first step
l.init <- (T*log(k))^0.5
lambda.1.cv <- 0.001*seq(l.init, 9*l.init, l.init/2)
blocks <- seq(0, T, floor(T/(floor(1*sqrt(T)))))
n.new <- length(blocks) - 1
blocks.size <- sapply(c(1:n.new), function(jjj) blocks[jjj+1] - blocks[jjj] )
bbb <- floor(n.new/5)
aaa <- sample(1:5, 1)
mu.sel <- 2                     ### we use fixed tuning parameter to estimate low-rank component

initial.phi <- matrix(0, k, k*p*(n.new))       ### this is the inital value for S_j's in our L+S_j paper which you can set from previous observations
cv.index <- seq(aaa, n.new, floor(n.new/bbb))  ### cv index, for cross validation index.

### estimating low-rank components
fit.initL <- fista.nuclear(X.sel, Y.sel, lambda = mu.sel, k, niter = 25, backtracking = TRUE, 0*diag(k))
L.new <- fit.initL$phi.hat
print(paste("estimated rank", qr(L.new)$rank))

#######################################################################################
### removing estimated low-rank effects from original data
dat.sel <- matrix(0, T, k)
dat.sel[1,] <- data[1,]
X.tmp <- data[1:(T-1),]
Y.tmp <- data[2:T,]
dat.sel[2:T,] <- Y.tmp - X.tmp %*% t(L.new)

### Apply sparse fused-lasso to estimate change point candidates by using the removed low-rank data
temp.first <- first.step.cv.new.blocks.rank("LASSO", dat.sel, weight = NULL, lambda.1.cv, 1, max.iteration = max.iteration,
                                            tol = tol, step.size = (1/2)*10^(-4), cv.index, blocks = blocks, initial = initial.phi)
print(temp.first$cv)
fisrt.brk.points <- temp.first$brk.points
pts <- fisrt.brk.points
print(paste("estimated change points:", pts))

pts.final <- pts
m.hat <- length(pts.final)
pts.est <- rep(0, m.hat+2)
pts.est[1] <- 0
pts.est[m.hat+2] <- T
pts.est[2:(m.hat+1)] <- pts.final

### Partition dataset according to the estimated change points
data.seg <- vector("list", m.hat+1)
for(j in 1:(m.hat+1)){
    data.seg[[j]] <- as.matrix(data[((pts.est[j]+1):pts.est[j+1]),])
}

### Setting parameters for estimating transition matrices at each segments
L.hat <- L.new
S.hat <- vector("list", m.hat+1)
for(j in 1:(m.hat+1)){ 
    S.hat[[j]] <- 0*diag(k)
}

loss.val <- 0
loss.val.new <- 1e8
max.iter = 10
rr <- 1
lambda.sel <- rep(0, m.hat+1)
obj.vals <- c()

lambda.vec <- vector("list", m.hat+1)
lambda.vec[[1]] <- seq(0.01, 0.05, length.out = 5)
lambda.vec[[2]] <- seq(0.01, 0.05, length.out = 5)
lambda.vec[[3]] <- seq(0.01, 0.05, length.out = 5)
lambda.vec[[4]] <- seq(0.01, 0.05, length.out = 5)
lambda.vec[[5]] <- seq(0.01, 0.05, length.out = 5)
lambda.vec[[6]] <- seq(0.01, 0.05, length.out = 5)
lambda.vec[[7]] <- seq(0.01, 0.05, length.out = 5)
lambda.vec[[8]] <- seq(0.01, 0.05, length.out = 5)
lambda.vec[[9]] <- seq(0.01, 0.05, length.out = 5)
lambda.vec[[10]] <- seq(0.01, 0.05, length.out = 5)
lambda.vec[[11]] <- seq(0.01, 0.05, length.out = 5)
lambda.vec[[12]] <- seq(0.01, 0.05, length.out = 5)

### Use initial low-rank component to initialize sparse components for each segments
for(j in 1:(m.hat+1)){
    len <- dim(data.seg[[j]])[1]
    X.seg <- data.seg[[j]][1:(len-1),]
    Y.seg <- data.seg[[j]][2:len,]
    fit.sparse <- bic.sparse(X.seg, Y.seg, lambda.vec[[j]], L.hat)
    S.hat[[j]] <- fit.sparse$S.hat
    print(paste("Estimating segment:", j, "...... with tuning", fit.sparse$lambda))
}

while(abs(loss.val.new - loss.val) > tol){
    if(rr > max.iter){
        break
    }else{
        print(paste(as.integer(rr), loss.val.new))
        rr <- rr + 1
        loss.val <- loss.val.new
        obj.vals <- c(obj.vals, loss.val)
        ### Use estimated (initialized) sparse component to re-estimate low-rank component
        ### Re-build dataset by removing sparse components
        Y.seg <- vector("list", m.hat+1)
        data.new.seg <- vector("list", m.hat+1)
        for(j in 1:(m.hat+1)){
            data.new.seg[[j]] = Y.seg[[j]] <- matrix(0, dim(data.seg[[j]])[1], dim(data.seg[[j]])[2])
            data.new.seg[[j]][1,] <- data.seg[[j]][1,]
            Y.seg[[j]] <- data.seg[[j]][-1,]
            for(jj in 1:(dim(data.seg[[j]])[1]-1)){
                xtm <- matrix(data.new.seg[[j]][jj,], 1, k)
                ytm <- matrix(Y.seg[[j]][jj,], 1, k)
                tmp <- ytm - xtm %*% t(S.hat[[j]])
                data.new.seg[[j]][jj+1,] <- tmp
            }
        }
        data.sparse <- c()
        for(j in 1:(m.hat+1)){
            data.sparse <- rbind(data.sparse, as.matrix(data.new.seg[[j]]))
        }
        X.lr <- data.sparse[1:(T-1),]
        Y.lr <- data.sparse[2:T,]
        fit.L <- fista.nuclear(X.lr, Y.lr, lambda = mu.sel, k, niter = 25, backtracking = TRUE, 0*diag(k))
        L.hat <- t(fit.L$phi.hat)
        print(qr(L.hat)$rank)
        
        ### Re-estimate sparse components until converge
        dat.tmp <- matrix(0, T, k); dat.tmp[1,] <- pre.data[1,]
        Y.tmp <- matrix(0, T-1, k); Y.tmp <- pre.data[(2:T),]
        for(j in 1:(T-1)){
            xtm <- matrix(dat.tmp[j,], 1, k)
            ytm <- matrix(Y.tmp[j,], 1, k)
            tmp <- ytm - xtm %*% t(L.hat)
            dat.tmp[j+1,] <- tmp
        }
        
        data.est <- vector("list", m.hat+1)
        X.sp = Y.sp <- vector("list", m.hat+1)
        for(j in 1:(m.hat+1)){
            data.est[[j]] <- as.matrix(dat.tmp[((pts.est[j]+1):pts.est[j+1]),])
            len <- dim(data.est[[j]])[1]
            X.sp[[j]] <- data.est[[j]][1:(len-1),]
            Y.sp[[j]] <- data.est[[j]][2:len,]
            fit.sparse <- bic.sparse(X.sp[[j]], Y.sp[[j]], lambda.vec[[j]], L.hat)
            S.hat[[j]] <- fit.sparse$S.hat
            lambda.sel[j] <- fit.sparse$lambda
            print(paste("Estimating segment:", j, "...... with tuning", fit.sparse$lambda))
        }
        
        ### stopping criterion values
        loss.val.new <- obj.func(L.hat, S.hat, pre.data[1:(T-1),], pre.data[2:T,], lambda = lambda.sel[1], mu = mu.sel)
        loss.val.new <- round(loss.val.new, 3)
        print("===============================================")
    }
}

### Put all estimated transition matrices together
phi.est <- matrix(0, k, k*(m.hat+1))
for(j in 1:(m.hat+1)){
    phi.est[,((j-1)*k+1):(k*j)] <- as.matrix(L.hat + S.hat[[j]])
}

### Plotting all components, low-rank and all sparse components
print(plot.matrix(abs(L.hat),1, name = "Low-rank"))
print(qr(L.hat)$rank)
for(j in 1:(m.hat+1)){
    print(plot.matrix(abs(S.hat[[j]]), 1, name = paste(j, "sparse")))
}

### Plotting the structure of sparse components for each segments
library(igraph)
for(seg in 1:(m.hat+1)){
    adj.mat <- matrix(0, k, k)
    for(i in 1:k){
        for(j in 1:k){
            if(S.hat[[seg]][i,j] != 0){
                adj.mat[i,j] <- sign(S.hat[[seg]][i,j])
                #adj.mat[i,j] <- S.hat[[seg]][i,j]
            }
        }
    }
    colnames(adj.mat) <- colnames(pre.data)
    net <- graph.adjacency(adj.mat, "directed", diag = FALSE)
    print(paste(seg, "edges:", length(E(net)) / (k^2)))
    # l <- layout_in_circle(net)
    # plot.igraph(net, vertex.label = V(net)$name, layout = l, 
    #             vertex.label.color = "black", edge.color = "black", edge.arrow.size = 0.1,
    #             vertex.shape ="none", vertex.color = "orange", vertex.label.cex = 0.7)
}

########### Plotting the jumps heatmap and corresponding pixels ###########
p <- dim(S.hat[[1]])[1]
thres <- rep(0, 12)
for(j in 1:12){
    idx <- order(S.hat[[j]], decreasing = TRUE)[1:20]
    thres[j] <- min(S.hat[[j]][idx])
}

# thres <- c(0.02, 0.0075, 0.01, 0.01, 0.0075, 0.0075, 0.01, 0.015, 0.02, 0.01, 0.0075, 0.0075)
par(mar = c(1, 1, 1, 1))
for(j in 1:12){
    bs <- matrix(0, nrow = 24, ncol = 32)
    pixels <- which(S.hat[[j]] >= thres[j], arr.ind = TRUE)
    for(i in 1:dim(pixels)[1]){
        pixels_ind <- pixels[i,]
        from <- c(pixels_ind[1] %% 24, (pixels_ind[1] %/% 24)+1)
        to <- c(pixels_ind[2] %% 24, (pixels_ind[2] %/% 24)+1)
        bs[from[1], from[2]] <- S.hat[[j]][pixels_ind[1], pixels_ind[2]]
        bs[to[1], to[2]] <- S.hat[[j]][pixels_ind[1], pixels_ind[2]]
    }
    x <- rescale(bs, c(-0.5, 0.5))
    #cellcol <- color.scale(cbind(x, c(-1, rep(1, 31))), c(0, 1), 0, c(1, 0))[,1:32]
    color2D.matplot(x, cs1 = c(1,0), cs2 = c(1,0), cs3 = c(1,0), xlab = "", ylab = "", axes = FALSE)
}