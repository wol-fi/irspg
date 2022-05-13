irspg <- function(X0, objfn=NULL, objfn_g=NULL, projfn=NULL, eqfn=NULL, eqfn_g=NULL, nmls=FALSE,
                    cntrl=list(eps_opt=1e-2, eta_min=1e-15, eta_max=1e+15, 
                    eps_feas=1e-6, kmax=50, ktmax=5, gtol=1e-01, ftol=1e-01), ...){
  gtol <- cntrl$gtol
  ftol <- cntrl$ftol
  eps_opt <- cntrl$eps_opt
  eta_min <- cntrl$eta_min
  eta_max <- cntrl$eta_max
  eps_feas <- cntrl$eps_feas
  eq_error <- eps_feas
  eps_opt <- cntrl$eps_opt
  kmax <- cntrl$kmax
  ktmax <- cntrl$ktmax
  
  # functions 
  etak_f <- function(){
    sk <- xk - xk_prev
    uk <- c(objfn_g(xk)) - c(objfn_g(xk_prev))
    
    su <- c(t(sk) %*% uk)
    if(su > 0){
      etak <- min(eta_max, max(eta_min, c(t(sk) %*% sk)/su ))
    } else{
      etak <- eta_max
    }
    return(etak)
  }
  
  restoration_fn <- function(xk){
    yhat <- projfn(xk)
    yk_f <- function(lambda) lambda*xk + (1-lambda)*yhat
    l <- 0.1
    ii <- 1 - seq(0, 1, length.out=100)
    for(ll in ii){
      yk <- yk_f(ll)
      if(abs(eqfn(yk)) <= rr*nC_xk){ 
        break
      }
    }
    return(yk)
  }
  
  # Initial Parameters
  
  xk <- projfn(X0)
  th_prev <- 0.99 # in (0,1)
  rr <- 0.25 # in [0,1)
  beta <- n*10 # > 0
  
  # Iterations
  xk_prev <- xk - 0.1
  
  nC_xk <- abs(eqfn(xk))
  objfn0 <- objfn(xk)
  
  ct <- Sys.time()
  for(k in 0:kmax){
    
    gr_prev <- norm(objfn_g(xk), "2")
    objfn_prev <- objfn(xk)
    
    # Step 1: Restoration, yk
    yk <- restoration_fn(xk)
    
    # Step 2
    omega_k <- 1/(k+1)*0.5
    thki_prev <- min(1, min(th_prev) + omega_k)
    
    # Step 3-4: Non-Monotone Line Search 
    etak <- etak_f()
    
    if(nmls){
      # non-monotone line search:
      for(kt in 1:ktmax){
        z <- yk + dk
        if(abs(eqfn(z)) <= rr*nC_xk){ #& norm(z - xk, "2") <= beta*nC_xk){
          yk <- z
          dk <- projfn(yk - etak * objfn_g(yk)) - yk
        } else {
          break
        }
      }
    } else {
      dk <- projfn(yk - etak * objfn_g(yk)) - yk  
    }

    nC_yk <- abs(eqfn(yk)) 
    
    # Step 5
    t <- 1
    gr_fyk <- objfn_g(yk)
    cdf <- crossprod(dk, gr_fyk)
    while(objfn(yk + t*dk) > (objfn(yk) + 0.1 * t * cdf)){
      t <- 0.9 * t
    }
    tdec <- tki <- t
    zki <- yk + tki*dk
    
    # alternative stopping criteria
    # if((abs(eqfn(zki)) <= eps_feas) & abs(gr_prev - sum(abs(objfn_g(zki)))) <= gtol){
    if(abs(eqfn(zki)) <= eps_feas){
      if(gr_prev -  norm(objfn_g(zki), "2") <= gtol | (objfn_prev - objfn(zki)) <= ftol){
        xk <- zki
        ct <- Sys.time() - ct
        return(list(X=xk, iter_out=k+1, eq_error=eqfn(xk), ctime=ct, objfn=objfn(xk), opti_at_first=(objfn(xk > objfn0))))
      }
    }
    
    # Step 6
    th <- seq(0, thki_prev, length.out = 5)
    dfz <- objfn(xk) - objfn(zki)
    dcy <- nC_xk - norm(eqfn(yk), "2")
    dcz <- nC_xk - norm(eqfn(zki), "2")
    
    aa <- th*dfz + (1-th)*dcy
    aa <- aa[aa >= 0.5*dcy]
    thki <- aa[which.max(aa)]
    
    xk_prev <- xk
    xk <- zki
    nC_xk <- abs(eqfn(xk))
    thk <- thki
    th_prev <- c(th_prev, thk)
  }
}








