# load("data/test_data.RData")
# source("R/example_fun.R")
# 
# 
# b <- sig %*% w
# vmb <- vm - t(b) %*% b
# fk <- 1
# x0 <- f_x0_2(A)
# B <- (b %*% t(b)) * J
# vex <- (vm - t(b) %*% b)
# 
# mdl <- irspgm(x0=x0, fn=fn, grad_fn=grad_fn, Ceq=Ceq, b=b, vmb=vmb, fk=1)
# 
# res$frob[i] <- fn(mdl$xk)
# res$time[i] <- mdl$ctime
# res$iter_out[i] <- mdl$iter_out
# res$vtol[i] <- mdl$vtol