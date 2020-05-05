####################################################################################################
################ This script demonstrates an example for a real applications #######################
####################################################################################################
rm(list = ls(all = TRUE))
data <- read.csv("Weekly.csv", header = TRUE)

####### Pre-processing data by removing NA columns and set up time stamps #######
library(dplyr)
time.stamp <- as.POSIXlt(data[,1], format = "%m/%d/%Y")
pre.data <- as.matrix(data[,-1] %>% 
                          select_if(~ !any(is.na(.))))

####################################################################################################
############################ Using algorithm to analyze change points ##############################
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

######## First step: Fused lasso change point detection, find candidates #########
### Variable related parameters
T <- dim(pre.data)[1]
k <- dim(pre.data)[2]
X.sel <- as.matrix(pre.data[1:(T-1),])
Y.sel <- as.matrix(pre.data[2:T,])

### Generating clusters to parallel computing
cl <- makeCluster(2)
registerDoParallel(cl)
getDoParWorkers()

####################################
######## Setting parameters ########
####################################
p = p.t <- 1                ### true AR lag
tol.sel <- 0.1
rep <- 0
tol <- 5*10^(1)             ### tolerance 
step.size <- 2*10^(-4)      ### step size 
max.iteration <- 100        ### max number of iteration for the Lasso solution
max.iter <- 10
method <- c('LASSO')

### In this application, we use fixed tuning parameter for low-rank component
# mu.sel <- 2.5     ### 1-factor
# mu.sel <- 0.75    ### 3-factor
mu.sel <- 0.50    ### 5-factor

### Tuning parameter for first step
l.init <- (T*log(k))^0.5
lambda.1.cv <- 0.0001*seq(l.init, 9*l.init, l.init/2)
blocks <- seq(0, T, floor(T/(floor(2*sqrt(T)))))
n.new <- length(blocks) - 1
blocks.size <- sapply(c(1:n.new), function(jjj) blocks[jjj+1] - blocks[jjj] )
bbb <- floor(n.new/5)
aaa <- sample(1:5, 1)

initial.phi <- matrix(0, k, k*p*(n.new))       ### this is the inital value for S_j's in our L+S_j paper which you can set from previous observations
cv.index <- seq(aaa, n.new, floor(n.new/bbb))  ### cv index, for cross validation index.

### Initialization for low-rank component
fit.initL <- fista.nuclear(X.sel, Y.sel, lambda = mu.sel, k, niter = 100, backtracking = TRUE, 0*diag(k))
L.new <- fit.initL$phi.hat
qr(L.new)$rank
print(plot.matrix(abs(L.new), 1))

###### Step 1 ######
### removing estimated low-rank effects from original data
dat.sel <- matrix(0, T, k)
dat.sel[1,] <- pre.data[1,]
X.tmp <- pre.data[1:(T-1),]
Y.tmp <- pre.data[2:T,]
dat.sel[2:T,] <- Y.tmp - X.tmp %*% t(L.new)

### Apply sparse fused-lasso to estimate change point candidates by using the removed low-rank data
temp.first <- first.step.cv.new.blocks.rank("LASSO", dat.sel, weight = NULL, lambda.1.cv, 1, max.iteration = max.iteration, 
                                            tol = tol, step.size = (1/2)*10^(-4), cv.index, blocks = blocks, initial = initial.phi)
print(temp.first$cv)
fisrt.brk.points <- temp.first$brk.points
pts <- fisrt.brk.points
print(fisrt.brk.points)

###### Step 2 ######
dat.screen <- matrix(0, T, k); dat.screen[1,] <- pre.data[1,]
X.tmp <- pre.data[1:(T-1),]
Y.tmp <- pre.data[2:T,]
dat.screen[-1,] <- Y.tmp - X.tmp %*% t(L.new)

n <- T - 1
## 1/40 for 1-factor
## 1/20 for 3-factor
## 1/15 for 5-factor
omega <- (1/20)*(((log(n))^1)*log(k))  ### the penalty term in the information criterion: the higher omega, the smaller number of break points selected.
lambda.2 <- (1/2)*(log(n)*log(k))/n    ### the second tuning parameter.
temp <- second.step(dat.screen, lambda = lambda.2, p, max.iteration = 200, tol = tol, step.size = 10^(-3), pts, omega)
final.brk.points <- temp$pts           ### final break points selected

MTSplot(pre.data)
abline(v = final.brk.points)

###### Step 3 ######
### Estimated change points
pts.final <- final.brk.points
# selected change points: 56, 84, 112, 336, 364, 392, 420, 504, 672, 756
m.hat <- length(pts.final)
pts.est <- rep(0, m.hat+2)
pts.est[1] <- 0
pts.est[m.hat+2] <- T
pts.est[2:(m.hat+1)] <- pts.final

### Partition dataset according to the estimated change points
data.seg <- vector("list", m.hat+1)
for(j in 1:(m.hat+1)){
    data.seg[[j]] <- as.matrix(pre.data[((pts.est[j]+1):pts.est[j+1]),])
}

###### Initialize parameters
L.hat <- L.new
S.hat <- vector("list", m.hat+1)
for(j in 1:(m.hat+1)){ 
    S.hat[[j]] <- 0*diag(k)
}

###### Starting Algorithm according to Theorem 4. (simplified version) ######
# loss.val <- 0
# loss.val.new <- 1e5
# rr <- 1
# lambda.sel <- rep(0, m.hat+1)
# obj.vals <- c()

### Tuning parameters
#(old version low-rank)
# mu.sel <- 0.15 ### for 1-factor model
# mu.sel <- 0.355 ### for 3-factor model.
# mu.sel <- 0.5 ### for 5-factor model.


## sparse tunings 
lambda.vec <- vector("list", m.hat+1)

### for 1-factor:
# lambda.vec[[1]] <- seq(0.0075, 0.015, length.out = 5)
# lambda.vec[[2]] <- seq(0.005, 0.01, length.out = 5)
# lambda.vec[[3]] <- seq(0.0025, 0.005, length.out = 5)

### for 3-factor: (full-algorithm version)
# lambda.vec[[1]] <- seq(0.0075, 0.01, length.out = 5)
# lambda.vec[[2]] <- seq(0.005, 0.01, length.out = 5)
# lambda.vec[[3]] <- seq(0.0025, 0.005, length.out = 5)
# lambda.vec[[4]] <- seq(0.005, 0.01, length.out = 5)
# lambda.vec[[5]] <- seq(0.005, 0.0075, length.out = 5)
# lambda.vec[[6]] <- seq(0.0075, 0.025, length.out = 5)
# lambda.vec[[7]] <- seq(0.01, 0.025, length.out = 5)
# lambda.vec[[8]] <- seq(0.005, 0.0075, length.out = 5)
# lambda.vec[[9]] <- seq(0.005, 0.0075, length.out = 5)
# lambda.vec[[10]] <- seq(0.0025, 0.005, length.out = 5)
# lambda.vec[[11]] <- seq(0.0025, 0.005, length.out = 5)

### for 3-factor:
lambda.vec[[1]] <- seq(0.01, 0.025, length.out = 5)
lambda.vec[[2]] <- seq(0.01, 0.025, length.out = 5)
lambda.vec[[3]] <- seq(0.0025, 0.0075, length.out = 5)
lambda.vec[[4]] <- seq(0.01, 0.05, length.out = 5)
lambda.vec[[5]] <- seq(0.0025, 0.0075, length.out = 5)
lambda.vec[[6]] <- seq(0.0025, 0.0075, length.out = 5)
lambda.vec[[7]] <- seq(0.01, 0.025, length.out = 5)
lambda.vec[[8]] <- seq(0.01, 0.05, length.out = 5)
lambda.vec[[9]] <- seq(0.0075, 0.025, length.out = 5)
lambda.vec[[10]] <- seq(0.005, 0.0075, length.out = 5)
lambda.vec[[11]] <- seq(0.005, 0.01, length.out = 5)

### for 5-factor: 
lambda.vec[[1]] <- seq(0.01, 0.05, length.out = 5)
lambda.vec[[2]] <- seq(0.01, 0.05, length.out = 5)
lambda.vec[[3]] <- seq(0.0025, 0.0075, length.out = 5)
lambda.vec[[4]] <- seq(0.01, 0.05, length.out = 5)
lambda.vec[[5]] <- seq(0.0025, 0.0075, length.out = 5)
lambda.vec[[6]] <- seq(0.0025, 0.0075, length.out = 5)
lambda.vec[[7]] <- seq(0.01, 0.025, length.out = 5)
lambda.vec[[8]] <- seq(0.01, 0.05, length.out = 5)
lambda.vec[[9]] <- seq(0.005, 0.025, length.out = 5)
lambda.vec[[10]] <- seq(0.005, 0.01, length.out = 5)
lambda.vec[[11]] <- seq(0.005, 0.01, length.out = 5)

### Estimate transition matrices for every segments
for(j in 1:(m.hat+1)){
    len <- dim(data.seg[[j]])[1]
    X.seg <- data.seg[[j]][1:(len-1),]
    Y.seg <- data.seg[[j]][2:len,] - X.seg %*% t(L.hat)            #### use removed data
    fit.sparse <- bic.sparse(X.seg, Y.seg, lambda.vec[[j]], L.hat)
    S.hat[[j]] <- fit.sparse$S.hat
    print(paste("Estimating segment:", j, "...... with tuning", fit.sparse$lambda))
}
