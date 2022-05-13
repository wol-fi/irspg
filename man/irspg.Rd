\name{irspg}
\alias{irspg}
\title{Spectral Projected Gradient with Inexact Restoration}
\description{An algorithm for large-scale nonlinear optimization problems with nonconvex constraints, based on the Spectral Projected Gradient method (SPG) with Inexacte Restoration (IR). The basic SPG is highly efficient for solving very large problems when projections on the feasbile set (constraints) are easy. The IR addon extends the algorithm for general nonlinearly constrained problems.}

\usage{irspgm(X0, objfn=NULL, objfn_g=NULL, projfn=NULL, eqfn=NULL, eqfn_g=NULL,nmls=FALSE, cntrl, ...)}

\arguments{
  \item{X0}{a numeric vector or matrix of initial values for the parameters to be optimized over.}
  \item{objfn}{the objective function to be minimized. Should be of class \code{function} and return a scalar result.}
  \item{objfn_g}{gradient of the objective function. Should be of class \code{function}.}
  \item{projfn}{projection function that projects \code{X} on the feasible set.}
  \item{eqfn}{equality constraint. Should be of class \code{function}.}
  \item{eqfn_g}{gradient of the equality constraint. Should be of class \code{function}.}
  \item{nmls}{Boolean indicating if non-monotone line-search should be used. This can speed up convergence.}
  \item{cntrl}{Optional: parameters to control the convergence.}
  \item{...}{additional values passed to the objective function }
}
\details{}
\value{Returns a list containing the optimal values of \code{X}, the final objective \code{objfn}, the accuracy of the equality constraint \code{eq_error} and other information.}
\references{Gomes-Ruggiero, M. A., Martinez, J. M., & Santos, S. A. (2009). Spectral projected gradient method with inexact restoration for minimization with nonconvex constraints. SIAM Journal on Scientific Computing, 31(3), 1628-1652.}