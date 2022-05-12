

# New algo 2021-05-04 -----------------------------------------------------

# Note: 
# the IR handles the equality constraint separately from the Omega constraint
# Omega represents X^2 <= 1
# C(X) defines the Lagrangian
# -> for projection, use the one as in the un-market-constrained case

irspg <- function(X0, fn=NULL, grad_fn=NULL, Ceq=NULL, grad_Ceq=NULL,
                    cntrl=list(eps_opt=1e-2, eta_min=1e-15, eta_max=1e+15, 
                    eps_feas=1e-6, kmax=50, ktmax=5, gtol=1e-01, ftol=1e-01), ...){
  gtol <- cntrl$gtol
  ftol <- cntrl$ftol
  eps_opt <- cntrl$eps_opt
  eta_min <- cntrl$eta_min
  eta_max <- cntrl$eta_max
  eps_feas <- cntrl$eps_feas
  vtol <- eps_feas
  eps_opt <- cntrl$eps_opt
  kmax <- cntrl$kmax
  ktmax <- cntrl$ktmax
  
  # fk <- dim(as.matrix(X0))[2]
  
  ## functions ##
  etak_f <- function(){
    sk <- xk - xk_prev
    uk <- c(grad_fn(xk)) - c(grad_fn(xk_prev))
    
    su <- c(t(sk) %*% uk)
    if(su > 0){
      etak <- min(eta_max, max(eta_min, c(t(sk) %*% sk)/su ))
    } else{
      etak <- eta_max
    }
    return(etak)
  }
  
  restoration_fn <- function(xk){
    yhat <- proj_f(xk)
    yk_f <- function(lambda) lambda*xk + (1-lambda)*yhat
    l <- 0.1
    ii <- 1 - seq(0, 1, length.out=100)
    for(ll in ii){
      yk <- yk_f(ll)
      if(abs(Ceq(yk)) <= rr*nC_xk){ 
        break
      }
    }
    return(yk)
  }
  
  ## Initial Parameters ##
  
  xk <- proj_f(X0)
  th_prev <- 0.99 # in (0,1)
  rr <- 0.25 # in [0,1)
  beta <- n*10 # > 0
  
  ## Iterations ##
  xk_prev <- xk - 0.1
  
  nC_xk <- abs(Ceq(xk))
  fn0 <- fn(xk)
  
  ct <- Sys.time()
  for(k in 0:kmax){
    # print(paste("iteration:", k, "| Gradient:", round(norm(grad_fn(xk), "2"), 3), "| Constraint:", round(Ceq(xk), 6), "| Objective:", round(fn(xk), 3)))
    # print(paste("iteration:", k, "    |dk|:", norm(dk, "2")))
    
    gr_prev <- norm(grad_fn(xk), "2")
    fn_prev <- fn(xk)
    
    # Step 1
    # Note: for equality constraints this writes different than in the paper,
    # see: https://www.ime.unicamp.br/~martinez/irexpanded.pdf
    
    # Step 1: Restoration, yk
    yk <- restoration_fn(xk)
    
    # Step 2
    omega_k <- 1/(k+1)*0.5
    thki_prev <- min(1, min(th_prev) + omega_k)
    
    # Step 3-4: Non-Monotone Line Search 
    etak <- etak_f()
    dk <- proj_f(yk - etak * grad_fn(yk)) - yk
    
    # # non-monotone line search:
    # for(kt in 1:ktmax){
    #   z <- yk + dk
    #   if(abs(Ceq(z)) <= rr*nC_xk){ #& norm(z - xk, "2") <= beta*nC_xk){
    #     yk <- z
    #     dk <- proj_f(yk - etak * grad_fn(yk)) - yk
    #   } else {
    #     break
    #   }
    # }
    # 
    nC_yk <- abs(Ceq(yk)) 
    
    # Step 5
    t <- 1
    # t <- 10
    gr_fyk <- grad_fn(yk)
    cdf <- crossprod(dk, gr_fyk)
    while(fn(yk + t*dk) > (fn(yk) + 0.1 * t * cdf)){
      t <- 0.9 * t
    }
    tdec <- tki <- t
    zki <- yk + tki*dk
    
    # alternative stopping criteria
    # if((abs(Ceq(zki)) <= eps_feas) & abs(gr_prev - sum(abs(grad_fn(zki)))) <= gtol){
    if(abs(Ceq(zki)) <= eps_feas){
      if(gr_prev -  norm(grad_fn(zki), "2") <= gtol | (fn_prev - fn(zki)) <= ftol){
        xk <- zki
        ct <- Sys.time() - ct
        return(list(xk=xk, iter_out=k+1, vtol=Ceq(xk), ctime=ct, fn=fn(xk), opti_at_first=(fn(xk > fn0))))
      }
    }
    
    # # Step 6
    th <- seq(0, thki_prev, length.out = 5)
    dfz <- fn(xk) - fn(zki)
    dcy <- nC_xk - norm(Ceq(yk), "2")
    dcz <- nC_xk - norm(Ceq(zki), "2")
    
    aa <- th*dfz + (1-th)*dcy
    aa <- aa[aa >= 0.5*dcy]
    thki <- aa[which.max(aa)]
    
    xk_prev <- xk
    xk <- zki
    nC_xk <- abs(Ceq(xk))
    thk <- thki
    th_prev <- c(th_prev, thk)
  }
}








