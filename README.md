# Modified and Truncated Newton Methods for Large-Scale Unconstrained Optimisation

Team project for **Numerical Optimization for Large Scale Problems** – MSc in Mathematical Engineering, Politecnico di Torino (A.Y. 2023/2024).

**Authors:** Lucia Ghezzi, Emanuela Girlando, Giuliana Paradiso

**Full report:** [Report.pdf](Report.pdf)

## Overview

MATLAB implementation and comparison of two Newton-type methods for unconstrained optimisation, both with a backtracking line search based on the Armijo sufficient-decrease condition:

- **Modified Newton method:** the Hessian is corrected as B_k = ∇²f(x_k) + τ_k I, with τ_k increased until a Cholesky factorisation succeeds, so that B_k is positive definite.
- **Truncated Newton method:** the Newton system is solved inexactly with the conjugate gradient (CG) method, stopped on a relative-residual tolerance or when a direction of negative curvature appears. Implemented both **without and with preconditioning** (incomplete Cholesky factorisation, with a diagonal preconditioner as fallback).

## Test problems

| Problem | Dimension | Starting points |
|---|---|---|
| Rosenbrock function | n = 2 | 5 |
| Generalized Broyden tridiagonal function (Problem 5) | n = 1000 | 4 |
| Penalty function 1 (Problem 27) | n = 1000 | 5 |
| Variably dimensioned function (Problem 28) | n = 1000 | 5 |

Gradients and Hessians are implemented analytically. For the large-scale problems the Hessian structure is exploited, for example the pentadiagonal Hessian of Problem 5.

## Parameters

| Parameter | Value |
|---|---|
| Armijo constant c1 | 1e-4 |
| Backtracking reduction factor ρ | 0.8 (0.5 for preconditioned Truncated Newton) |
| Maximum backtracking steps | 10 |
| Gradient-norm tolerance | 1e-8 |
| Maximum outer iterations | 5000 |

## Evaluation

For each method, problem and starting point the code reports the number of iterations, the final gradient norm, the execution time and the experimental rate of convergence. The rate is computed from three consecutive errors, using the known minimiser for the Rosenbrock function and the difference between successive iterates otherwise.

## Main results

- **57 runs in total, 54 converged** to the required tolerance. The 3 failures (Problem 28) stopped at gradient norms of order 1e-7 to 1e-8.
- The **Modified Newton method** was generally the most robust and efficient.
- The **Truncated Newton method without preconditioning** was faster on some starting points of Problems 27 and 28.
- **Preconditioning** increased execution times in most cases and was decisive for convergence in only one case (Problem 28, first starting point).

Full tables and discussion are in the [report](Report.pdf).

## Repository structure

```
.
├── README.md
├── Report.pdf
└── UO_assignment_code.m   # main script and all method implementations
```

## How to run

Open `UO_assignment_code.m` in MATLAB and run it section by section. Each test problem has its own block of sections (Modified Newton, Truncated Newton without and with preconditioning). Change the starting point by uncommenting the corresponding `x0` line at the top of each block.

## Notes

The code is published as submitted for the exam, so it is consistent with the results in the report.
