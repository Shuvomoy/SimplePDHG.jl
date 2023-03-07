using LinearAlgebra

# Write a julia struct that will take c which is a vector of length n, A which is a matrix of size m x n, and b which is a vector of length m

struct LP_Data{T<:AbstractFloat, I <: Integer}
    c::Vector{T} # cost vector of length n
    A::Matrix{T} # date matrix of size m x n
    b::Vector{T} # resource vector of length m
    m::I # number of rows of A
    n::I # number of columns of A
end

# contains the problem setting, i.e., the parameter γ as set by the user

struct PDHG_settings

    # user settings to solve the problem using PDHG 
    # =============================================
    η::Float64 # first step size of PDHG
    τ::Float64 # second step size of PDHG
    maxit::Int64 # maximum number of iteration
    tol::Float64 # tolerance, i.e., if |||| ≤ tol, we take x to be an optimal solution
    verbose::Bool # whether to print information about the iterates
    freq::Int64 # how often print information about the iterates

    # constructor for the structure, so if user does not specify any particular values, 
    # then we create a setting object with default values
    function PDHG_settings(η, τ, maxit, tol, verbose, freq) 
        new(η, τ, maxit, tol, verbose, freq)
    end
    
end

# construct PDHG state that is going to contain the necessary information to completely describe one iteration of PDHG algorithm

mutable struct PDHG_state{T <: AbstractVecOrMat{<: Real}, I <: Integer} # contains information regarding one iterattion sequence
    x::T # iterate x_n
    y::T # iterate y_n
    η::Float64 # step size
    τ::Float64 # step size
    k::I # iteration counter  
end

## We need a method to construct the initial value of the type PDHG_state

function PDHG_state(problem::LP_Data)
    n = problem.n
    m = problem.m
    σmaxA = norm(problem.A,2)
    η_preli = (1/(σmaxA)) - 1e-6
    τ_preli = (1/(σmaxA)) - 1e-6
    x_0 = zeros(n)
    y_0 = zeros(m)
    return PDHG_state(x_0, y_0, η_preli, τ_preli, 1)
end

## Write one iteration of PDHG

function PDHG_iteration!(problem::LP_Data, state::PDHG_state)

    # unpack the current state information
    x_k = state.x
    y_k = state.y
    k = state.k

    # compute the next iterate

    # compute x_k_plus_1 = x_k - η*problem.A'*y_k- η*c
    x_k_plus_1 = x_k - state.η*(problem.A'*y_k)- state.η*problem.c

    # project on to the positive orthant
    project_nonnegative!(x_k_plus_1)

    # compute y_k_plus_1 = y + τ*A*(2*x_k_plus_1 - x_k) - τ*b
    y_k_plus_1 = y_k + state.τ*problem.A*(2*x_k_plus_1 - x_k) - state.τ*problem.b

    # load the computed values in the PDHG_state
    state.x = x_k_plus_1
    state.y = y_k_plus_1
    state.k = k + 1

    # return the updated state
    return state

end