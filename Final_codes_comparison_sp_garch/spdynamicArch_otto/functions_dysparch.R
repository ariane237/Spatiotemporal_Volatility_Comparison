GMM_SDPD_2SLS_ARCH_ind_timelags <- function(Y, X, W, info){
  
  ksy = info$ksy # spatial exapansion order of y
  ksx = info$ksx # spatial exapansion order of x
  
  if(is.null(Y)){
    stop("Y is missing")
  }
  # if(is.null(X)){
  # stop("X is missing")
  # }
  if(is.null(W)){
    stop("W is missing")
  }
  if(is.null(info)){
    stop("info is missing")
  }
  
  dimW <- dim(W)
  n    <- dimW[1]
  if(length(dimW) == 2){
    p <- 1
    new.W <- array(, dim = c(n,n,p))
    new.W[,,1] <- W
    W <- new.W
  } else {
    p <- dimW[3]
  }
  
  if(length(dimW) > 3 |  dimW[1] !=  dimW[2] | dimW[2] !=  n){
    stop("Check W matrix (W must be of dimension n x n x p)")
  }
  
  if(dim(Y)[1] != n | length(dim(Y)) != 2){
    stop("Y should be a matrix of dimension n x T")
  }
  
  t <- dim(Y)[2] - 1
  
  # cat("Number of cross-sectional units:", n,  "\n")
  # cat("Length of the time series:", t,  "\n")
  # cat("Number of spatial lags (weight matrices):", p,  "\n")
  
  
  s  <- t - 1
  nt <- n * t
  ns <- n * s
  
  yt   <- as.vector(Y[, 2:(t+1)])
  ytlv <- as.vector(Y[, 1:(t)]) # ytl vector
  ytl  <- array(0, dim = c(nt,n)) # ytl matrix to get individual temporal coefficients
  ysl  <- array(0, dim = c(nt,p))
  ystl <- array(0, dim = c(nt,p))
  
  for(i in 1:n){
    ytl[seq(1, nt, by = n) + i - 1, i] <- Y[i, 1:(t)]
  }
  
  for(i in 1:t){
    for(j in 1:p){
      ysl[(1+(i-1)*n):(i*n), j]  <- W[,,j] %*% yt[(1+(i-1)*n):(i*n)];
      ystl[(1+(i-1)*n):(i*n), j] <- W[,,j] %*% ytlv[(1+(i-1)*n):(i*n)];
    }
  }
  
  if(info$stl + info$tl == 2){
    xw <- cbind(ytl, ystl);
  } else if (info$stl + info$tl == 1){
    if(info$stl == 1){
      xw <- ystl
    } else {
      xw <- ytl
    }
  } else if (info$stl + info$tl == 0){
    stop("No spatial and no temporal lag given")
  } else {
    stop("Double-Check stl & tl # in Info structure")
  }
  
  X_stacked <- NULL
  W1hx <- NULL
  MHX <- NULL
  if(!is.null(X)){
    for(i in 1:t){
      X_stacked <- rbind(X_stacked, X[,i+1,])
    }
    
    if(dim(X)[3] == 1){
      X_stacked <- array(X_stacked, dim = c(length(X_stacked), 1))
    }
  }
  
  xs  <- X_stacked;
  xt  <- cbind(xw, xs);
  zt  <- cbind(ysl, xt);
  
  kz  <- dim(zt)[2]
  kx  <- dim(xt)[2]
  kxs <- dim(xs)[2]
  kxw <- dim(xw)[2]
  
  c <- sqrt((t-(1:s))/(t-(1:s)+1)); 
  
  F <- diag(t)
  F <- F[, 1:(t-1)]
  for(i in 1:(t-1)){
    F[(i+1):t, i] <- -1/(t-i);
    F[, i]        <- c[i] * F[,i];
  }
  
  hyt  <- array(yt, dim = c(n,t));
  hyt  <- hyt %*% F;
  hyt  <- as.vector(hyt);
  hytlv <- array(ytlv, dim = c(n,t));
  hytlv <- hytlv %*% F;
  hytlv <- as.vector(hytlv);
  
  hytl  <- array(0, dim = c(ns,n)) # hytl matrix to get individual temporal coefficients
  for(i in 1:n){
    hytl[seq(1, ns, by = n) + i - 1, i] <- array(hytlv, dim = c(n,s))[i, ]
  }
  
  hysl <- array(ysl, dim = c(n,t,p));
  hysltemp <- array(0, dim = c(n,t-1,p));
  for(i in 1:p){
    hysltemp[,,i] <- hysl[,,i] %*% F;
  }
  hysl <- array(hysltemp, dim = c(ns,p));
  
  hystl <- array(ystl, dim = c(n,t,p));
  hystltemp <- array(0, dim = c(n,t-1,p));
  for(i in 1:p){
    hystltemp[,,i] <- hystl[,,i] %*% F;
  }
  hystl <- array(hystltemp, dim = c(ns,p));
  
  if(!is.null(X_stacked)){
    kx <- dim(X_stacked)[2]
    hx <- array(X_stacked, dim = c(n,t,kx));
    hxtemp <- array(0, dim = c(n,t-1,kx));
    for(i in 1:kx){
      hxtemp[,,i] <- hx[,,i] %*% F;
    }
    hx <- array(hxtemp, dim = c(ns,kx));
  } else {
    hx <- NULL
  }
  
  if(info$stl + info$tl == 2){
    hxw <- cbind(hytl, hystl);
  } else if(info$stl + info$tl == 1){
    if(info$stl == 1){
      hxw <- hystl
    } else {
      hxw <- hytl
    } 
  } else if(info$stl + info$tl == 0){
    stop("no spatial and temporal lag given, check suitability")
  } else {
    stop("Doube-check info$stl and info$tl")
  }
  
  pyid     <- array(0, dim = c(ksy,2));
  pa       <- 1
  pb       <- p
  pyid[1,] <- c(pa, pb);
  for(k in 2:ksy){
    pa <- pa + p^(k-1)
    pb <- pb + p^k
    pyid[k,] <- c(pa, pb);
  }
  
  WY <- array(0, dim = c(nt, pyid[ksy,2]));
  WY[, 1:p] = ystl
  for(i in 1:t){
    for(k in 1:(ksy-1)){
      for(j in 1:p){
        WY[(1+(i-1)*n):(i*n), (pyid[k,2] + 1 + (j-1)*p^k):(pyid[k,2]+j*p^k)] <- W[,,j] %*% WY[(1+(i-1)*n):(i*n), pyid[k,1]:pyid[k,2]];
      }
    }
  }
  
  if(!is.null(X_stacked)){
    kx <- dim(X_stacked)[2]
    W1hx <- array(0, dim = c(ns,kx*p))
    
    for(i in 1:s){
      for(j in 1:p){
        W1hx[(1+(i-1)*n):(i*n), (1+(j-1)*kx):(j*kx)] <- W[,,j] %*% hx[(1+(i-1)*n):(i*n),];
      }
    }
    
    pxid <- array(0, dim = c(ksx,2))
    pa       <- 1
    pb       <- p*kx
    pxid[1,] <- c(pa, pb);
    for(k in 2:ksx){
      pa <- pa + p^(k-1)*kx
      pb <- pb + p^k*kx
      pxid[k,] <- c(pa, pb);
    }
    
    WHX <- array(0, dim = c(ns, pxid[ksx,2]*kx))
    WHX[,1:(p*kx)] <- W1hx;
    for(i in 1:s){
      for(k in 1:(ksx-1)){
        for(j in 1:p){
          WHX[(1+(i-1)*n):(i*n), (pxid[k,2]+1+(j-1)*p^k*kx):(pxid[k,2]+j*p^k*kx)] <- W[,,j] %*% WHX[(1+(i-1)*n):(i*n), pxid[k,1]:pxid[k,2]];
        }
      }
    }
  }
  
  
  ## Following is the IV without interaction term
  
  pyid <- array(0, dim = c(ksy,2));
  pa   <- 1;
  pb   <- p;
  pyid[1,] <- c(pa, pb);
  for(k in 2:ksy){
    pa <- pa + p
    pb <- pb + p
    pyid[k,] <- c(pa, pb)
  }
  
  MY <- array(0, dim = c(nt, pyid[ksy,2]))
  MY[,1:p] <- ystl;
  
  for(i in 1:t){
    for(k in 1:(ksy-1)){
      for(j in 1:p){
        MY[(1+(i-1)*n):(i*n), (pyid[k,2]+j):(pyid[k,2]+j)] <- W[,,j] %*% MY[(1+(i-1)*n):(i*n), (pyid[k,1]-1+j):(pyid[k,1]-1+j)];
      }
    }
  }
  
  if(!is.null(X_stacked)){
    
    kx <- dim(X_stacked)[2]
    pxid <- array(0, dim = c(ksx,2))
    pa <- 1
    pb <- p*kx
    pxid[1,] <- c(pa, pb)
    for(k in 2:ksx){
      pa <- pa + p*kx
      pb <- pb + p*kx
      pxid[k, ] <- c(pa, pb)
    }
    
    MHX <- array(0, dim = c(ns,pyid[ksx,2]*kx));
    MHX[,1:(p*kx)] <- W1hx;
    for(i in 1:s){
      for(k in 1:(ksx-1)){
        for(j in 1:p){
          MHX[(1+(i-1)*n):(i*n),(pxid[k,2]+1+(j-1)*kx):(pxid[k,2]+j*kx)] <- W[,,j] %*% MHX[(1+(i-1)*n):(i*n),(pxid[k,1]+(j-1)*kx):(pxid[k,1]-1+j*kx)];
        }
      }
    }
    
  }
  
  Qw     <- cbind(ytl, WY);
  Qw_alt <- cbind(ytl, MY);
  
  if(ksx == 0){
    Qs     <- NULL
    Qs_alt <- NULL
    qs     <- 0
    qs_alt <- 0
  } else {
    Qs     <- cbind(hx, W1hx)
    Qs_alt <- cbind(hx, MHX)
    qs     <- dim(Qs)[2]
    qs_slt <- dim(Qs_alt)[2]
  }
  
  Qw     <- Qw[1:ns,]
  Qw_alt <- Qw_alt[1:ns,]
  qw     <- dim(Qw)[2]
  qw_alt <- dim(Qw_alt)[2]
  
  Q      <- cbind(Qw, Qs)
  Q_alt  <- cbind(Qw_alt, Qs_alt)
  kq     <- dim(Q)[2]
  kq_alt <- dim(Q_alt)[2]
  
  hxs <- hx
  hxt <- cbind(hxw, hxs)
  hzt <- cbind(hysl, hxt)
  
  Qhz <- array(0, dim = c(kq, kz))
  QQ  <- array(0, dim = c(kq, kq))
  Qhy <- array(0, dim = c(kq, 1))
  
  Qhz_alt <- array(0, dim = c(kq_alt, kz))
  QQ_alt  <- array(0, dim = c(kq_alt, kq_alt))
  Qhy_alt <- array(0, dim = c(kq_alt, 1))
  
  Jn <- diag(n) - array(1/n, dim = c(n,n));
  
  for(i in 1:s){
    if(info$ted == 1){
      
      Qhz <- Qhz + t(Q[(1+(i-1)*n):(i*n),]) %*% Jn %*% hzt[(1+(i-1)*n):(i*n),];
      QQ  <- QQ  + t(Q[(1+(i-1)*n):(i*n),]) %*% Jn %*% Q[(1+(i-1)*n):(i*n),];
      Qhy <- Qhy + t(Q[(1+(i-1)*n):(i*n),]) %*% Jn %*% hyt[(1+(i-1)*n):(i*n)];
      
      Qhz_alt <- Qhz_alt + t(Q_alt[(1+(i-1)*n):(i*n),]) %*% Jn %*% hzt[(1+(i-1)*n):(i*n),];
      QQ_alt  <- QQ_alt  + t(Q_alt[(1+(i-1)*n):(i*n),]) %*% Jn %*% Q_alt[(1+(i-1)*n):(i*n),];
      Qhy_alt <- Qhy_alt + t(Q_alt[(1+(i-1)*n):(i*n),]) %*% Jn %*% hyt[(1+(i-1)*n):(i*n)];
      
    } else {
      
      Qhz <- Qhz + t(Q[(1+(i-1)*n):(i*n),]) %*% hzt[(1+(i-1)*n):(i*n),];
      QQ  <- QQ  + t(Q[(1+(i-1)*n):(i*n),]) %*% Q[(1+(i-1)*n):(i*n),];
      Qhy <- Qhy + t(Q[(1+(i-1)*n):(i*n),]) %*% hyt[(1+(i-1)*n):(i*n)];
      
      Qhz_alt <- Qhz_alt + t(Q_alt[(1+(i-1)*n):(i*n),]) %*% hzt[(1+(i-1)*n):(i*n),];
      QQ_alt  <- QQ_alt  + t(Q_alt[(1+(i-1)*n):(i*n),]) %*% Q_alt[(1+(i-1)*n):(i*n),];
      Qhy_alt <- Qhy_alt + t(Q_alt[(1+(i-1)*n):(i*n),]) %*% hyt[(1+(i-1)*n):(i*n)];
      
    }
  }
  
  theta <- mldivide(t(Qhz) %*% ginv(QQ) %*% Qhz, t(Qhz) %*% ginv(QQ) %*% Qhy);
  theta_alt <- mldivide(t(Qhz_alt) %*% ginv(QQ_alt) %*% Qhz_alt, t(Qhz_alt) %*% ginv(QQ_alt) %*% Qhy_alt);
  
  e <- hyt - hzt %*% theta;
  e_alt <- hyt - hzt %*% theta_alt;
  if(info$ted == 1){
    for(i in 1:s){
      e[(1+(i-1)*n):(i*n)] <- Jn %*% e[(1+(i-1)*n):(i*n)];
      e_alt[(1+(i-1)*n):(i*n)] <- Jn %*% e_alt[(1+(i-1)*n):(i*n)];
    }
  }
  
  sigma2 <- mean((e-mean(e))^2);
  sigma4 <- mean((e-mean(e))^4);
  
  sigma2_alt <- mean((e_alt-mean(e_alt))^2);
  sigma4_alt <- mean((e_alt-mean(e_alt))^4);
  
  lambda <- theta[1:p];
  delta  <- theta[(p+1):kz];
  
  lambda_alt <- theta_alt[1:p];
  delta_alt  <- theta_alt[p+1:kz];
  
  lambdaW     <- array(0, dim = c(n,n));
  lambdaW_alt <- array(0, dim = c(n,n));
  for(j in 1:p){
    lambdaW     <- lambdaW     + lambda[j] * W[,,j];
    lambdaW_alt <- lambdaW_alt + lambda_alt[j] * W[,,j];
  }
  Sn <- diag(n) - lambdaW;
  Sn_alt <- diag(n) - lambdaW_alt;
  
  DSiD <- 1/(sigma2*ns) * t(Qhz) %*% ginv(QQ) %*% Qhz;
  DSiD_alt <- 1/(sigma2_alt*ns) * t(Qhz_alt) %*% ginv(QQ_alt) %*% Qhz_alt;
  
  SIG <- tryCatch(1/ns * solve(DSiD), error = function(e){cat("DSiD not invertible \n"); return(array(NA, dim = dim(DSiD)))})
  # SIG <- 1/ns * solve(DSiD);
  
  std <- sqrt(abs(diag(SIG)));
  tstat <- theta/std;
  
  SIG_alt <- tryCatch(1/ns * solve(DSiD_alt), error = function(e){cat("DSiD_alt not invertible \n"); return(array(NA, dim = dim(DSiD_alt)))})
  # SIG_alt <- 1/ns * solve(DSiD_alt);
  std_alt <- sqrt(abs(diag(SIG_alt)));
  tstat_alt <- theta_alt/std_alt;
  
  results <- list(theta = theta, std = std, SIG = SIG, tstat = tstat, sigma2 = sigma2,
                  theta_alt = theta_alt, std_alt = std_alt, SIG_alt = SIG_alt, tstat_alt = tstat_alt, sigma2_alt = sigma2_alt,
                  e = array(e, dim = c(n, t)), e_alt = array(e_alt, dim = c(n, t)), hyt = array(hyt, dim = c(n, t)))
  
  return(results)
  
}

mldivide <- function(A, b){
  return(ginv(t(A) %*% A) %*% t(A) %*% b)
  # return(qr.solve(A, b))
}
