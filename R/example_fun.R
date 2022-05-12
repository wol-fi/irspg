# f_x0_2 <- function(A){
#   # modified version
#   e <- eigen(A)
#   vv <- e$vectors[,1:fk,drop=F]
#   if(mean(vv) < 0) vv <- -vv
#   ll <- e$values[1:fk]
#   aopt <- c()
#   for(j in 1:fk){
#     lam <- ll[j]
#     v <- vv[,j]
#     a1 <- sqrt((lam-1) * norm(v, "2")^2 / (norm(v, "2")^4 - sum(v^4)))
#     a2 <- 1/(sqrt(max(abs(v))))
#     aopt <- c(aopt, c(min(a1, a2)))
#   }
#   if(fk == 1){
#     x0 <- vv * c(aopt)
#   } else {
#     x0 <- vv %*% diag(aopt)  
#   }
#   
#   sel <- which(colMeans(x0)<0)
#   x0[,sel] <- -x0[,sel]
#   return(x0)
# }
# 
# # equality constraint
# Ceq <- function(x){
#   x <- matrix(x, n)
#   c(t(w) %*% sig %*% (x %*% t(x) * J) %*% sig %*% w - vmb)
# }
# 
# # gradient of equality constraint
# gr_C <- function(x){
#   x <- matrix(x, n)
#   c(- 2* sig %*% (J * w %*% t(w)) %*% sig %*% x)
# }
# 
# # objective function
# fn <- function(x, Ahat){
#   x <- matrix(x, n)
#   norm(J * x %*% t(x)  - Ahat, "F")^2  
# }
# 
# # gradient of objective fun
# grad_fn <- function(x){
#   x <- matrix(x, n)
#   c(4*(J * x %*% t(x) - Ahat) %*% x)
# }
# 
# # projection function, general
# proj_f <- function(y){
#   y <- matrix(y, n)
#   for(kk in 1:100){
#     y <- Peq(y)
#     rs <- rowSums(y^2) 
#     if(max(rs) > 1.02){
#       ind <- sqrt(rs) 
#       ind[ind < 1] <- 1
#       ind <- ind + (ind - 1) * 0.1
#       y <- y/ind
#     } else {
#       break
#     }
#     if(abs(Ceq(y)) < vtol) break
#   }
#   return(c(y))
# }
# 
# # projection function, equality constraint
# Peq <- function(Y){
#   XX <- Y%*%t(Y)
#   aa <- c(t(b) %*% (B %*% XX %*% B * J) %*% b)
#   bb <- 2*c(t(b)%*%((XX %*% B)*J) %*% b)
#   cc <- c(t(b) %*% (XX * J) %*% b - vex)
#   dd <- max(0, bb^2 - 4*aa*cc)
#   lam <- (-bb + sqrt(dd))/(2*aa)
#   return((I + lam*B) %*% Y)
# }