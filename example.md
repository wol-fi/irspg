## The Task

Find the nearest correlation matrix subject to an equality constrained.
This problem was described in Schadner and Traut (2022) at the
estimation of forward-looking stock correlation matrices.

## Solving the Problem

The test sample covers implied volatilities, weights, Pearson’s
correlations and pre-computed values for the cross-section and index of
the S&P 100 as observed on January 31, 1996.

#### 1. Defining the Functions


    ## Objective Function
    # i.e., the function to be minimized
    objfn <- function(x){
      x <- matrix(x, n)
      norm(J * x %*% t(x)  - Ahat, "F")^2
    }

    # gradient of objective function
    objfn_g <- function(x){
      x <- matrix(x, n)
      c(4*(J * x %*% t(x) - Ahat) %*% x)
    }

    ## Equality Constraint
    # i.e., the C(x) function within Gomes-Ruggiero et al. (2009)
    eqfn <- function(x){
      x <- matrix(x, n)
      c(t(w) %*% sig %*% (x %*% t(x) * J) %*% sig %*% w - vmb)
    }

    # gradient of equality constraint
    eqfn_g <- function(x){
      x <- matrix(x, n)
      c(- 2* sig %*% (J * w %*% t(w)) %*% sig %*% x)
    }

    ## Projection Function
    # i.e., the projection on the tangible set
    # see Gomes-Ruggiero et al. (2009) for greater details
    projfn <- function(y){
      Peq <- function(Y){
        XX <- Y%*%t(Y)
        aa <- c(t(b) %*% (B %*% XX %*% B * J) %*% b)
        bb <- 2*c(t(b)%*%((XX %*% B)*J) %*% b)
        cc <- c(t(b) %*% (XX * J) %*% b - vex)
        dd <- max(0, bb^2 - 4*aa*cc)
        lam <- (-bb + sqrt(dd))/(2*aa)
        return((I + lam*B) %*% Y)
      }
      
      y <- matrix(y, n)
      for(kk in 1:100){
        y <- Peq(y)
        rs <- rowSums(y^2)
        if(max(rs) > 1.02){
          ind <- sqrt(rs)
          ind[ind < 1] <- 1
          ind <- ind + (ind - 1) * 0.1
          y <- y/ind
        } else {
          break
        }
        if(abs(eqfn(y)) < vtol) break
      }
      return(c(y))
    }

#### 2. Load the Data and Run the Optimization

    library(irspg)

    ## loading and assigning the sample data
    data(sp100) 
    tmp <- sapply(ls(), function(z) assign(z, get(z))) 

    ## Run the optimization
    mdl <- irspg(x0, objfn, objfn_g, projfn, eqfn, eqfn_g)


    ## Results

    # (preview) the best set of parameters found
    head(mdl$X) 
    #> [1] 0.7757936 0.5698422 0.6091928 0.6844730 0.5415583 0.6312939

    # final value of objective function
    mdl$objfn
    #> [1] 109.7039

    # starting value of objective function
    mdl$opti_at_first
    #> [1] 287.9409

    # accuracy of meeting the equality constraint
    mdl$eq_error
    #> [1] 5.20417e-18

    # computational time
    mdl$ctime
    #> Time difference of 0.05481315 secs

## References

Schadner, W & Traut, J. (2022). “Estimating Forward-Looking Stock
Correlations from Risk Factors”, *Mathematics* 10(10), 1649.
<https://doi.org/10.3390/math10101649>.
