## Bayesian multilabel modelling functions

## function to generate from multivariate normal
my.mvrnorm = function(n, mu, Sigma){
  p = length(mu)
  # compute square root of covariance matrix
  eo=eigen(Sigma, symmetric=TRUE)
  sigma.sqrt=eo$vec%*%diag(eo$val^0.5)%*%t(eo$vec)
  
  # generate random normals from runif by box-muller transform
  rnorm.vec = sqrt(-2*log(runif(n*p)))*cos(2*pi*runif(n*p))
  
  # generate sample matrix
  sample.matrix = matrix(rep(mu, n), nrow=n, byrow=T) +
    matrix(rnorm.vec, nrow=n, ncol=p)%*%sigma.sqrt
  return(sample.matrix)
}

etPois = function(lam) lam*exp(lam) / (exp(lam)-1)

rtpois <- function(n, lambda,tol=1e-10){
  ## Simulate from zero-truncated Poisson
  
  ## Initialize output
  x <- rep(NA,n)
  
  ## Identify lambda values below tolerance
  low <- which(lambda < tol)
  nlow <- length(low)
  
  if(nlow > 0){
    x[low] <- 1
    
    if(nlow < n)
      x[-low] <- qpois(runif(n-nlow, dpois(0, lambda[-low]), 1), lambda[-low])
  }
  else
    x <- qpois(runif(n-nlow, dpois(0, lambda), 1), lambda)
  
  return(x)
}

#### Using Gibbs sampler
bmlml.Gibbs = function(Y1, X1, K, nsamp=1e3, maxit=1e3, tol=1e-3){
  require(MCMCpack)
  require(BayesLogit)
  Y = t(Y1); X = t(X1)
  L = nrow(Y)
  N = ncol(Y)
  D = nrow(X)
  which.plus = which(Y==1, arr.ind=T)
  
  ## initialize quantities
  W.list = list()
  V.list = list()
  U.list = list()
  M.list = list()
  W.list[[1]] = matrix(.1, nrow=D, ncol=K)
  V.list[[1]] = matrix(.1, nrow=L, ncol=K)
  U.list[[1]] = matrix(.1, nrow=K, ncol=N)
  M.list[[1]] = array(.1, dim=c(L, K, N))
  M23 = apply(M.list[[1]], 2:3, sum)
  
  
  for(i in 2:nsamp){
    W = matrix(1, nrow=D, ncol=K)
    V = matrix(1, nrow=L, ncol=K)
    U = matrix(1, nrow=K, ncol=N)

    WtX = t(W.list[[i-1]]) %*% X
    P = 1 / (1 + exp( - WtX))
    # r = rgamma(K, shape = 0.01, rate = 0.01)
    r = rep(1, K)
    
    ## sample V and U
    M12 = apply(M.list[[i-1]], 1:2, sum)

    for(k in 1:K){
      V[,k] = rdirichlet(1, 1+M12[,k])
      U[k,] = rgamma(rep(1,N), r[k]+M23[k,], P[k,])
    }
    
    ## sample W
    M = array(0, dim=c(L, K, N))
    for(ind in 1:nrow(which.plus)){
      l.plus = as.numeric(which.plus[ind,1])
      n.plus = as.numeric(which.plus[ind,2])
      for(k in 1:K){
        M[l.plus, k, n.plus] = rtpois(1, V[l.plus, k] * U[k, n.plus])
      }
    }
    M23 = apply(M, 2:3, sum)
    
    for(k in 1:K){
      Omega.k = rpg.gamma(N, r[k] + M23[k,], WtX[k,])
      Sig.k = ginv(X %*% diag(Omega.k) %*% t(X) + diag(1,D))
      mu.k = Sig.k %*% X %*% (M23[k,] - r[k])/2
      W[,k] = my.mvrnorm(1, mu.k, Sig.k)
    }
    
    ## update
    W.list[[i]] = W
    V.list[[i]] = V
    U.list[[i]] = U
    M.list[[i]] = M
  }
  
}

#### Using EM algorithm
bmlml.EM = function(Y1, X1, K, maxit=1e3, tol=1e-3){
  require(MASS)
  
  L = ncol(Y1)
  N = nrow(Y1)
  D = ncol(X1)
  one <- rep(1, N)
  meanx <- drop(one %*% X1)/N
  xc <- scale(X1, meanx, FALSE)         # first subtracts mean
  normx <- sqrt(drop(one %*% (xc^2)))
  names(normx) <- NULL
  xs <- scale(xc, FALSE, normx)        # now rescales with norm (not sd)
  
  Y = t(Y1); X = t(xs)
  which.plus = which(Y==1, arr.ind=T)
  
  ## initialize quantities
  W = matrix(1, nrow=D, ncol=K)
  V = matrix(1, nrow=L, ncol=K)
  eU = matrix(1, nrow=K, ncol=N)
  M = matrix(1, nrow=K, ncol=N)

  ## initialize hyperparameters
  r = rgamma(K, shape = 0.01, rate = 0.01)
  tau = rgamma(D, 0.01, 0.01)
  
  ## iterate
  eps.vec = rep(0 ,maxit)
  for(iter in 1:maxit){
    # r = rgamma(K, shape = 0.01, rate = 0.01)
    r1 = matrix(r, nrow=K, ncol=N, byrow=T)
    WtX = t(W) %*% X
    P = 1 / (1 + exp( - WtX))
    
    ## E step
    eM.array = array(0, dim=c(L, K, N))
    for(ind in 1:nrow(which.plus)){
      l.plus = as.numeric(which.plus[ind,1])
      n.plus = as.numeric(which.plus[ind,2])
      for(k in 1:K){
        eM.array[l.plus, k, n.plus] = V[l.plus, k] * eU[k, n.plus]
      }
    }
    eM.array = etPois(eM.array)
    eM.array[which(is.na(eM.array), arr.ind=T)] = 0
    eM.array[which(eM.array==Inf, arr.ind=T)] = 0
    eM23 = apply(eM.array, 2:3, sum)
    
    eU = (1/P) * (eM23 + r1)
    eOmega = ((eM23 + r1) / (2*WtX)) * tanh(WtX/2)
    
    ## M step
    V1 = matrix(0, nrow=L, ncol=K)
    W1 = matrix(0, nrow=D, ncol=K)
    eM12 = apply(eM.array, 1:2, sum)
    for(k in 1:K){
      V1[,k] = eM12[,k] / sum(eM12[,k])
      # tau = rgamma(D, 0.01, 0.01)
      Sk = X %*% diag(eOmega[k,]) %*% t(X) + diag(tau)
      dk = X %*% (eM23[k,] - r[k])/2
      W1[,k] = ginv(Sk) %*% dk
    }
    
    ## check convergence
    eps = norm(V1-V)/norm(V) + norm(W1-W)/norm(W)
    eps.vec[iter] = eps
    if(eps <= tol){
      break
    }
    else{
      V = V1; W = W1
    }
  }
  
  converged = TRUE
  if(iter==maxit & eps > tol){
    converged = FALSE
  }
  
  model = list(V=V, W=W, r=r, eU=eU, converged=converged, eps.vec=eps.vec[1:iter])
  class(model) = "BMLML"
  model
}

bmlml.EM1 = function(Y1, X1, K, r.tau, maxit=1e3, tol=1e-3){
  require(MASS)
  
  L = ncol(Y1)
  N = nrow(Y1)
  D = ncol(X1)

  Y = t(Y1); X = t(X1)
  which.plus = which(Y==1, arr.ind=T)
  
  ## initialize quantities
  W = matrix(1, nrow=D, ncol=K)
  V = matrix(1, nrow=L, ncol=K)
  eU = matrix(1, nrow=K, ncol=N)
  M = matrix(1, nrow=K, ncol=N)
  
  ## initialize hyperparameters
  r = r.tau[1:K]
  tau = r.tau[-(1:K)]
  # tau = rgamma(D,1,1)

  ## iterate
  eps.vec = rep(0 ,maxit)
  for(iter in 1:maxit){
    # r = rgamma(K, shape = 0.01, rate = 0.01)
    r1 = matrix(r, nrow=K, ncol=N, byrow=T)
    WtX = t(W) %*% X
    P = 1 / (1 + exp( - WtX))
    
    ## E step
    eM.array = array(0, dim=c(L, K, N))
    for(ind in 1:nrow(which.plus)){
      l.plus = as.numeric(which.plus[ind,1])
      n.plus = as.numeric(which.plus[ind,2])
      for(k in 1:K){
        eM.array[l.plus, k, n.plus] = V[l.plus, k] * eU[k, n.plus]
      }
    }
    eM.array = etPois(eM.array)
    eM.array[which(is.na(eM.array), arr.ind=T)] = 0
    eM.array[which(eM.array==Inf, arr.ind=T)] = 0
    eM23 = apply(eM.array, 2:3, sum)
    
    eU = (1/P) * (eM23 + r1)
    eOmega = ((eM23 + r1) / (2*WtX)) * tanh(WtX/2)
    
    ## M step
    V1 = matrix(0, nrow=L, ncol=K)
    W1 = matrix(0, nrow=D, ncol=K)
    eM12 = apply(eM.array, 1:2, sum)
    for(k in 1:K){
      V1[,k] = eM12[,k] / sum(eM12[,k])
      # tau = rgamma(D, 0.01, 0.01)
      Sk = X %*% diag(eOmega[k,]) %*% t(X) + diag(tau)
      dk = X %*% (eM23[k,] - r[k])/2
      W1[,k] = ginv(Sk) %*% dk
    }
    
    ## check convergence
    eps = norm(V1-V)/norm(V) + norm(W1-W)/norm(W)
    eps.vec[iter] = eps
    if(eps <= tol){
      break
    }
    else{
      V = V1; W = W1
    }
  }
  
  converged = TRUE
  if(iter==maxit & eps > tol){
    converged = FALSE
  }
  
  model = list(V=V, W=W, r=r, eU=eU, converged=converged, eps.vec=eps.vec[1:iter])
  class(model) = "BMLML"
  model
}

bmlml.EM2 = function(Y1, X1, Z1, K, maxit=1e3, tol=1e-3){
  require(MASS)
  
  L = ncol(Y1)
  N = nrow(Y1)
  D = ncol(X1)
  E = ncol(Z1)
  one <- rep(1, N)
  meanx <- drop(one %*% X1)/N
  xc <- scale(X1, meanx, FALSE)         # first subtracts mean
  normx <- sqrt(drop(one %*% (xc^2)))
  names(normx) <- NULL
  xs <- scale(xc, FALSE, normx)        # now rescales with norm (not sd)
  
  Y = t(Y1); X = t(xs); Z = t(Z1)
  which.plus = which(Y==1, arr.ind=T)
  
  ## initialize quantities
  W = matrix(1, nrow=D, ncol=K)
  S = matrix(1, nrow=E, ncol=K)
  V = matrix(1, nrow=L, ncol=K)
  eU = matrix(1, nrow=K, ncol=N)
  M = matrix(1, nrow=K, ncol=N)
  
  ## initialize hyperparameters
  r = rgamma(K, shape = 0.01, rate = 0.01)
  s = rgamma(K, shape = 0.01, rate = 0.01)
  tau = rgamma(D, 0.01, 0.01)
  eta = rgamma(E, 0.01, 0.01)
  
  ## iterate
  eps.vec = rep(0 ,maxit)
  for(iter in 1:maxit){
    # r = rgamma(K, shape = 0.01, rate = 0.01)
    r1 = matrix(r, nrow=K, ncol=N, byrow=T)
    s1 = t(matrix(s, nrow=K, ncol=L, byrow=T))
    WtX = t(W) %*% X
    P = 1 / (1 + exp( - WtX))
    StZ = t(Z) %*% S
    Q = 1 / (1 + exp( - StZ))
    
    ## E step
    eM.array = array(0, dim=c(L, K, N))
    for(ind in 1:nrow(which.plus)){
      l.plus = as.numeric(which.plus[ind,1])
      n.plus = as.numeric(which.plus[ind,2])
      for(k in 1:K){
        eM.array[l.plus, k, n.plus] = V[l.plus, k] * eU[k, n.plus]
      }
    }
    eM.array = etPois(eM.array)
    eM.array[which(is.na(eM.array), arr.ind=T)] = 0
    eM.array[which(eM.array==Inf, arr.ind=T)] = 0
    eM23 = apply(eM.array, 2:3, sum)
    eM12 = apply(eM.array, 1:2, sum)
    
    eU = (1/P) * (eM23 + r1)
    eOmega = ((eM23 + r1) / (2*WtX)) * tanh(WtX/2)
    eZeta = ((eM12 + s1) / (2*StZ)) * tanh(StZ/2)

    ## M step
    V1 = matrix(0, nrow=L, ncol=K)
    W1 = matrix(0, nrow=D, ncol=K)
    S1 = matrix(0, nrow=E, ncol=K)
    for(k in 1:K){
      V1[,k] = (eM12[,k] + Q[,k]) / sum(eM12[,k] + Q[,k])
      
      # tau = rgamma(D, 0.01, 0.01)
      Sk = X %*% diag(eOmega[k,]) %*% t(X) + diag(tau)
      dk = X %*% (eM23[k,] - r[k])/2
      W1[,k] = ginv(Sk) %*% dk
      
      Tk = Z %*% diag(eZeta[,k]) %*% t(Z) + diag(eta)
      ek = Z %*% (eM12[,k] - s[k])/2
      S1[,k] = ginv(Tk) %*% ek
    }
    
    ## check convergence
    eps = norm(V1-V)/norm(V) + norm(W1-W)/norm(W) + norm(S1-S)/norm(S)
    eps.vec[iter] = eps
    if(eps <= tol){
      break
    }
    else{
      V = V1; W = W1; S = S1
    }
  }
  
  converged = TRUE
  if(iter==maxit & eps > tol){
    converged = FALSE
  }
  
  model = list(V=V, W=W, S=S, r=r, eU=eU, converged=converged, eps.vec=eps.vec[1:iter])
  class(model) = "BMLML"
  model
}

predict.bmlml = function(model, newX, transform="scale", standardize=NULL){
  N = nrow(newX)
  L = nrow(model$V)
  
  if(!is.null(standardize)){
    newX = scale(newX, center=standardize[1,], scale=standardize[2,])
  } # scales the predictor variable if necessary

  ## get outputs
  Ypred = matrix(0, N, L)
  for(n in 1:N){
    for(l in 1:L){
      Ypred[n,l] = with(model, 1 - prod((V[l,] * exp(t(W) %*% newX[n,]) + 1)^(-r)))
    }
  }
  
  # return, applying transformation on predicted probabilities if necessary
  if(transform=="rank"){
    for(l in 1:L){
      Ypred[,l] = rank(Ypred[,l])/N
    }
  }else if(transform=="scale"){
    for(l in 1:L){
      Ypred[,l] = Ypred[,l] / max(Ypred[,l])
    }
  }
  Ypred
}

partial.bmlml = function(model, data, Y.index, X.index, lims=NULL, bins=1e2, median=F, ...){
  prob.mat = matrix(0, nrow=nrow(data), ncol=bins)
  
  ## set limits if not provided
  if(is.null(lims)){
    lims = c(min(data[,X.index]), max(data[,X.index]))
  }
  
  # get a sequence of values, replace each value in dataset and predict
  valueSeq = seq(lims[1], lims[2], length.out=bins)
  for(i in 1:bins){
    idata = data
    idata[,X.index] = valueSeq[i]
    ipred = predict.bmlml(model, idata, ...)
    prob.mat[,i] = ipred[,Y.index]
  }
  
  if(median){
    summary = apply(prob.mat,2,median)
  } else{
    summary = apply(prob.mat,2,mean)
  }
  cbind(valueSeq, summary)
}

selectParams = function(Y1, X1, K, nU, ...){
  L = ncol(Y1)
  N = nrow(Y1)
  D = ncol(X1)

  set.seed(10212016)
  r = runif(K+D,0,.2)
  model = bmlml.EM1(Y1, scale(X1), K=K, r.tau=r, ...)
  P= t(1/(1+exp(-scale(X1) %*% model$W)))
  P.fac = as.numeric(P/(1-P))
  P.fac[which(P.fac==Inf)] = 1e5 # control against dividing by 0
  
  # multiple sets of posterior draws... rows of U.mat
  U.mat = matrix(0, nrow=nU, ncol=N)
  for(num in 1:nU){
    m = rnegbin(rep(1,N),r[1],P)
    U = rgamma(rep(1,N),r[1]+m,P)
    U.mat[num,] = U
  }
  
  tuning.mat = matrix(0, nrow=1e3, ncol=K+D+1)
  for(i in 1:1e2){
    r1 = runif(K+D,0,.2)
    ratio.mat = dgamma(U.mat, r[1], P.fac)/dgamma(U.mat, r1[1], P.fac)
    bf = apply(ratio.mat, 2, function(x) prod(x, na.rm=T))
    
    tuning.mat[i,-1] = r1
    tuning.mat[i,1] = mean(bf[which(bf!=Inf)])
  }
  
  r.hat = tuning.mat[which.max(tuning.mat[,1]),-1]  
  r.hat
}