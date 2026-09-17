NEWS/Changelog

# quadrupen 1.1-0	(2026-09-17)

- performance
  - active set: no reallocation of the Gram matrices when variables enter or leave the model, LAPACK triangular solves, in-place Cholesky downdates
  - several variables (or groups) activated at once in the working set algorithm
    (new control option `maxadd`, default 10 for sparse models and 5 for group models)
  - proximal solvers: adaptive restart of FISTA, safer Lipschitz constant (Lanczos), step-invariant stopping rule and safeguarded Anderson acceleration for PGD
  - exact block coordinate descent with Newton steps for the `quadra` method of group models
  - ridge and lava: no p x p algebra for a diagonal structuring matrix (the default)
  - faster information criteria; cross-validation folds built on demand (less memory)
  - typical speed-ups over 1.0-0 with default settings: x2 to x4 for lasso, MCP and SCAD, x5 to x15 for group models, x9 to x36 for elastic-net and sparse designs with large active sets, x8 for ridge and x30 for lava with p = 8000; all lambdas now converge on these cases
- bug fixes (results may change)
  - robustness when the active set exceeds the rank of the design (no more failure of the Cholesky factorization or of the direct solve)
  - `bounded_reg()` with `method = "quadra"` (the default) stopped too early and could return suboptimal solutions; the optimality gap is now computed from the exact KKT conditions
  - l1/linf group penalty: wrong proximal operator when the group lies in the l1 ball (`fista` and `pgd` converged to wrong solutions)
  - `quadra` method for group models: exact solution of each block; sparse-group, coop and l1/linf variants now converge
  - proximal solvers could stop after one iteration on badly scaled data
  - convergence status now reflects the final optimality gap
  - the Lipschitz constant estimation no longer draws from R's random number generator
- minor updates
  - objects are rebuilt when a C++ header changes; `ARMA_NO_DEBUG` also on Windows

# quadrupen 1.0-0	(2026-06-05)

- major updates
  - complete rewriting of R code using R6 classes
  - complete rewriting of C++ code using template and OO style programming
  - included 'FusedLasso' from the archived package by Holger Hoefling (fixed CRAN's complaints)
  - added group-lasso/group-elastic and variants (group penalty: l1/l2, l1/linf, cooperative Lasso)
  - added sparse group-lasso/group-elastic and variant (group penalty: l1/l2, l1/linf, cooperative Lasso)
  - added lava and/post-lava (combination of sparse and dense regularization, Chernozukov et al, 2017)
  - extended to group-lava (group penalty: l1/l2, l1/linf, cooperative Lasso)
  - added mcp and scad (concave penalties, Chernozukov et al, 2017)
  - added refit version of Lasso/Elastic-Net ("relaxed" Lasso/Enet)
  - changing many parameters (badly) named, do not expect backward compatibility
  - added vignettes
- minor updates
  - Integration of changes from CRAN versions from 0.2-4 to 0.2-13
  - set up github workflow for pkgdown page
  - various fixes, more testing

# quadrupen 0.2-4	(2014-01-16)

  Minor:
    - memory leak corrected (sp_mat declaration)
    - linking to Rcpp/RcppArmadillo headers (requires R 3.0-2)

# quadrupen 0.2-4	(2013-11-01)
- added a 'lasso' function, simple wrapper to the elastic-net 'function'
- added computation of degrees of freedom (for elastic net and bounded regression)
- added a method to compute penalized criteria (BIC/AIC) of a quadrupen fit, with plot

# 0.2-3	(2013-08-26)

- added back the 'normalize' parameter
- standardization is performed within the C++ code
- use of sparse conversion from Matrix to Armadillo
- corrected bug with the 'intercept' and 'residuals' components of the quadrupen class
- added more tests in the inst directory
- correction in the documentation
- added r.squared to the quadrupen class

# quadrupen 0.2-2	(2013-04-08)

- minor fix to comply with recent ggplot2 updates.

# quadrupen 0.2-1	(2013-02-27)

- minor fix to pass CRAN check on Windows operating systems.

quadrupen 0.2-0	(2013-02-26)

- Major updates
  - added bounded regression (regression penalized by infinity norm + structured l2 norm)
  - added corresponding features for cross-validation and stability path
- Minor updates
- corrected wrong annotations of the stability path (PFER)
- handled normalization internally ('normalize' is no longer a parameter)
- more simple internal handling of penscales and correction of the rescaling of the intercept
- better use of multicore features
- handled runtime error exception in RcppArmadillo when the system is singular (end of the solution path)
  A consequence is quadrupen is less likely to crash due to user's "bad" parametrisation
- simplification of the C++ code, bugs corrected, probably new ones added :-'(
- added 'examples' and 'tests' directories

# quadrupen 0.1-0	(2012-10-09)

- first build: structured elastic-net with (weighted) quadratic loss, cross-validation and stability selection methods.
