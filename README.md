### *IR-SPG*
# Spectral Projected Gradient with Inexact Restoration
An algorithm for large-scale nonlinear optimization problems with nonconvex constraints, based on the Spectral Projected Gradient method (SPG) with Inexacte Restoration (IR). The basic SPG is highly efficient for solving very large problems when projections on the feasbile set (constraints) are easy. The IR addon extends the algorithm for general nonlinearly constrained problems.

#### Installation
``` r
    devtools::install_github("wol-fi/irspg")
```

#### Usage
``` r
irspg(X0, objfn, objfn_g, projfn, eqfn, eqfn_g)
```
where
* `X0` ... initial values
* `objfn` ... objective function
* `obfn_g` ... gradient of objective function
* `projfn` ... projection function
* `eqfn` ... equality constraint (function)
* `qfn_g` ... gradient of equality constraint

#### Example
Check the 'example.md' file for an illustrative instruction.

#### Keywords
non-linear programming, spectral projected gradient, inexact restoration, large-scale optimization

#### MSC Codes
90C30, 90C26, 65K05, 49K35

#### Applications
Estimation of Forward-Looking Stock Correaltion Matrices: <https://doi.org/10.3390/math10101649>

#### References
Gomes-Ruggiero, M. A., Martínez, J. M., & Santos, S. A. (2009). Spectral projected gradient method with inexact restoration for minimization with nonconvex constraints. *SIAM Journal on Scientific Computing, 31(3)*, 1628-1652. <https://doi.org/10.1137/070707828>
