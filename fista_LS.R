###### estimation functions for low-rank and sparse components

fista.nuclear <- function(A, b, lambda, d, niter, backtracking = TRUE, phi.true = diag(d)){
    #' This function is to estimate low-rank components via nuclear norm penalty. The main algorithm
    #' to solve the optimization problem is FISTA. 
    #' @param A: design matrix
    #' @param b: response vector 
    #' @param lambda: tuning parameter for low-rank component
    #' @param d: dimension of coefficient matrix
    #' @param niter: the number of iterations for FISTA to converge
    #' @param backtracking: A boolean argument, indicating whether employ backtracking in algorithm
    #' @param phi.true: true transition matrx, it's not necessary, default value is diagonal matrix
    #' @return A list including the estimated low-rank component, and objective function curves.
    #' @export
    tnew = t <- 1
    x <- matrix(0, d, d)
    xnew <- x
    y <- x
    AtA <- t(A) %*% A
    Atb <- t(A) %*% b
    
    obj.val = rel.err <- c()
    if(backtracking == TRUE){
        L <- norm(A, "2")^2 / 5
        eta <- 2
    }else{
        L <- norm(A, "2")^2
    }
    for(i in 1:niter){
        if(backtracking == TRUE){
            L.bar <- L
            flag <- FALSE
            while(flag == FALSE){
                prox <- prox.nuclear.func(y, A, b, L.bar, lambda, AtA, Atb)
                
                ### restricted part
                for(rr in 1:d){
                    for(cc in 1:d){
                        if(abs(prox[rr,cc]) > 0.25){
                            ### here we set alpha introduced in our paper equals to 0.25
                            prox[rr,cc] <- 0.25*sign(prox[rr,cc])
                        }
                    }
                }
                
                if(f.func(prox, A, b) <= Q.func(prox, y, A, b, L.bar, AtA, Atb)){
                    flag <- TRUE
                }else{
                    L.bar <- L.bar * eta
                }
            }
            L <- L.bar
        }
        x <- xnew
        xnew <- prox
        t <- tnew
        tnew <- (1 + sqrt(1 + 4*t^2)) / 2
        y <- xnew + ((t - 1) / tnew) * (xnew - x)
        
        obj.val <- c(obj.val, f.func(xnew, A, b) + nuclear.pen(xnew, lambda))
        rel.err <- c(rel.err, norm(xnew - phi.true, "F") / norm(phi.true, "F"))
    }
    return(list(phi.hat = xnew, obj.vals = obj.val, rel.err = rel.err))
}

fista.sparse <- function(A, b, lambda, d, niter, backtracking = TRUE, phi.true){
    #' This function is to compute the estimator of sparse component via L1 penalty. 
    #' @param A: design matrix
    #' @param b: response vector 
    #' @param lambda: tuning parameter for sparse component
    #' @param d: dimension of coefficient matrix
    #' @param niter: the number of iterations for FISTA to converge
    #' @param backtracking: A boolean argument, indicating whether employ backtracking in algorithm
    #' @param phi.true: true transition matrx, it's not necessary, default value is diagonal matrix
    #' @return A list including the estimated sparse component, and objective function curves.
    #' @export
    tnew = t <- 1
    x <- matrix(0, d, d)
    xnew <- x
    y <- x
    AtA <- t(A) %*% A
    Atb <- t(A) %*% b
    
    obj.val = rel.err <- c()
    if(backtracking == TRUE){
        L <- norm(A, "2")^2 / 5
        eta <- 2
    }else{
        L <- norm(A, "2")^2
    }
    for(i in 1:niter){
        if(backtracking == TRUE){
            L.bar <- L
            flag <- FALSE
            while(flag == FALSE){
                prox <- prox.sparse.func(y, A, b, L.bar, lambda, AtA, Atb)
                if(f.func(prox, A, b) <= Q.func(prox, y, A, b, L.bar, AtA, Atb)){
                    flag <- TRUE
                }else{
                    L.bar <- L.bar * eta
                }
            }
            L <- L.bar
        }
        x <- xnew
        xnew <- prox
        t <- tnew
        tnew <- (1 + sqrt(1 + 4*t^2)) / 2
        y <- xnew + ((t - 1) / tnew) * (xnew - x)
        
        obj.val <- c(obj.val, f.func(xnew, A, b) + sparse.pen(xnew, lambda))
        rel.err <- c(rel.err, norm(xnew - phi.true, "F") / norm(phi.true, "F"))
    }
    return(list(phi.hat = xnew, obj.vals = obj.val, rel.err = rel.err))
}

shrinkage <- function(y, tau){
    #' This function is to compute the soft-thresholding for lasso estimator
    #' @param y: input matrix
    #' @param tau: threshold value
    #' @export
    
    z <- matrix(0, nrow = nrow(y), ncol = ncol(y))
    for(i in 1:nrow(y)){
        for(j in 1:ncol(y)){
            z[i,j] <- sign(y[i,j]) * max(abs(y[i,j]) - tau, 0)
        }
    }
    return(z)
}

shrinkage.lr <- function(y, tau){
    #' This function is to compute the soft-thresholding for nuclear estimator: we threshold the singular
    #' value vector of input matrix
    #' @param y: input singular value vector of matrix. It's derived by SVD.
    #' @param tau: threshold value
    #' @return the thresholded singular value vector
    #' @export

    z <- rep(0, length(y))
    for(i in 1:length(y)){
        z[i] <- sign(y[i]) * max(0, abs(y[i]) - tau)
    }
    return(z)
}

f.func <- function(x, A, b){
    #' This function computes the quadratic loss
    #' @param x: coefficient matrix
    #' @param A: design matrix
    #' @param b: response vector
    #' @export
    
    return(0.5 * norm(A %*% x - b, "F")^2)
}

gradf.func <- function(x, AtA, Atb){
    #' This function is gradient function with respect to the quadratic loss
    #' @param x: coefficient matrix
    #' @param AtA: crossproduct of A
    #' @param Atb: crossproduct of A and b
    #' @export
    
    return(AtA %*% x - Atb)
}

nuclear.pen <- function(x, lambda){
    #' This function indicates the nuclear norm for low-rank penalty term
    #' @param x: coefficient
    #' @param lambda: tuning parameter
    #' @export
    
    d <- svd(x)$d
    return(lambda * sum(d))
}

sparse.pen <- function(x, lambda){
    #' This function indicates the L1 norm for sparse penalty term
    #' @param x: coefficient
    #' @param lambda: tuning parameter
    #' @export 

    return(lambda*sum(x))
}

Q.func <- function(x, y, A, b, L, AtA, Atb){
    #' This function is introduced by original FISTA paper, and it's used for updating the coefficient
    #' @param x: coefficient
    #' @param y: previous updated coefficient
    #' @param A: design matrix
    #' @param b: response matrix
    #' @param L: learning rate
    #' @param AtA: crossproduct of A
    #' @param Atb: crossproduct of A and b
    #' @export
    
    return(f.func(y, A, b) + sum((x - y) * gradf.func(y, AtA, Atb)) + 0.5 * L * norm(x - y, "F")^2)
}

prox.nuclear.func <- function(y, A, b, L, lambda, AtA, Atb){
    #' This function is to compute proximal function of low-rank component
    #' @param y: updated coefficient
    #' @param A: design matrix
    #' @param b: response matrix
    #' @param L: learning rate
    #' @param lambda: tuning parameter of low-rank component
    #' @param AtA: crossproduct of A
    #' @param Atb: crossproduct of A and b
    #' @export
    
    Y <- y - (1 / L) * gradf.func(y, AtA, Atb)
    d <- shrinkage.lr(svd(Y)$d, 2*lambda / L)
    return(svd(Y)$u %*% diag(d) %*% t(svd(Y)$v))
}

prox.sparse.func <- function(y, A, b, L, lambda, AtA, Atb){
    #' This function is to compute proximal function of sparse component
    #' @param y: updated coefficient
    #' @param A: design matrix
    #' @param b: response matrix
    #' @param L: learning rate
    #' @param lambda: tuning parameter of sparse component
    #' @param AtA: crossproduct of A
    #' @param Atb: crossproduct of A and b
    #' @export
    Y <- y - (1 / L) * gradf.func(y, AtA, Atb)
    return(shrinkage(Y, 2*lambda / L))
}

obj.func <- function(x.lr, x.sparse, A, b, lambda, mu){
    #' This function is the full objective function
    #' @param x.lr: low-rank coefficient
    #' @param x.sparse: sparse coefficients, in our paper, this is a list
    #' @param A: design matrix
    #' @param b: response matrix
    #' @param lambda: tuning parameter for sparse component
    #' @param mu: tuning parameter for low-rank component
    #' @export
    
    m <- length(x.sparse)
    loss <- 0
    for(i in 1:m){
        loss <- loss + f.func((x.lr + x.sparse[[i]]), A, b) + sparse.pen(x.sparse[[i]], lambda)
    }
    return(loss + nuclear.pen(x.lr, mu))
}
