# SimplePDHG.jl

I created this simple educational Julia package (less than 350 lines of code and less than 3 hours coding in total) to demonstrate to the students of the MIT Course 15.084/6.7220 Nonlinear Optimization that how easy it is to implement an algorithm in Julia and connecting it to the optimization modeling language `JuMP.jl` so that anyone can use your package. 

Big thanks to [Oscar Dowson](https://odow.github.io/) for providing `MathOptInterface.jl ` code to connect this simple solver to `JuMP`! ([discourse link](https://discourse.julialang.org/t/connecting-a-simple-first-order-solver-to-solve-standard-form-linear-program-to-jump/95694))

## What does `SimplePDHG.jl` do?

This is an educational package used to demonstrate the ease of implementing an algorithm in `Julia` and incorporating it with one of Julia's main optimization modeling language `JuMP`. The package  designed to solve linear programming problems of the form:

```julia
minimize    c'x
subject to  A x = b
            G h ≤ h
            x ∈ ℝ^n
```

where `x` is the decision variable. Under the hood the `SimplePDHG.jl` implements the vanilla PDHG algorithm (see Section 3.3 of [this book](https://large-scale-book.mathopt.com/LSCOMO.pdf)) to solve standard form linear optimization problem of the form `min{c'x ∣ Ax=b, x ≥ 0, x ∈ ℝ^n}`.

##  Installation 

```julia
] add https://github.com/Shuvomoy/SimplePDHG.jl.git
```

## Usage through `JuMP`

```julia
using JuMP
model =  Model(SimplePDHG.Optimizer)
@variable(model, x >= 0)
@variable(model, 0 <= y <= 3)
@objective(model, Min, 12x + 20y)
@constraint(model, c1, 6x + 8y >= 100)
@constraint(model, c2, 7x + 12y >= 120)
optimize!(model)
println("Objective value: ", objective_value(model))
println("x = ", value(x))
println("y = ", value(y))
```

Output should be:

```julia
Objective value: 205.000090068938
x = 14.999887019427522
y = 1.2500722917903861
```

## Vector syntax in JuMP

Thanks to `JuMP` and `MathOptInterface.jl `, we can use vectorized syntax to solve our optimization problem as well. 

```julia
# data 
A = [1 1 9 5; 3 5 0 8; 2 0 6 13]
b = [7, 3, 5]
c = [1, 3, 5, 2]
m, n = size(A)

# JuMP code
using JuMP
model =  Model(SimplePDHG.Optimizer)
@variable(model, x[1:n] >= 0)
@objective(model, Min, c'*x)
@constraint(model, A*x .== b)
optimize!(model)
println("Objective value: ", objective_value(model))
println("x = ", value.(x))
```

The output should be:

```julia
Objective value: 4.92307853474916

x = [0.4230758379471922, 0.3461545013289533, 0.6923078385630216, 0.0]
```



