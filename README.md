# SimplePDHG

Implements vanilla PDHG to solve standard form linear optimization problem:

min c'x
s.t. Ax = b
     x ≥ 0
     x ∈ ℝ^n.



## Installation 

```
] add https://github.com/Shuvomoy/SimplePDHG.jl.git
```

## Example usage

```
## Problem data

A = [0.24960819905908674 0.8788337376233512 0.43078688988461517; 0.06806539450164809 2.1454929611189417 -0.8066868870829861]

b = [0.5723911460577251, 0.523953352474369]

c = [0.7728142977622404, 0.860578976651018, 0.002363122195560588]

## Default Stepsize

σmaxA = norm(A,2)

η_preli = 1/(2*σmaxA^2)

τ_preli = (1/(2*σmaxA^2)) - 1e-6

m, n = size(A)

## Create LP data object

problem = LP_Data(c, A, b, m, n)

## Settings for the optimization solver

maxit_test = 1000000

tol_test = 1e-6

verbose_test = true

freq_test = 10000

settings = PDHG_settings(η_preli, τ_preli, maxit_test, tol_test, verbose_test, freq_test)

## Call the solver

final_state_PDHG = PDHG_solver(problem, settings)

## Optimal solution

x_star = final_state_PDHG.x
```



