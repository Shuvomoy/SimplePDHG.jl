# Given a vector x in R^n, write code to project on to the positive orthant {x ∣ x ≥ 0}. The projection should be done in place, i.e. the function should modify the input vector x. The function should return the number of negative entries in x before projection.

# project_nonnegative!(x)

function project_nonnegative!(x::AbstractVector{T}) where T<:Real
    for i in eachindex(x)
        x[i] = max(zero(T), x[i])
    end
end

# Write the same function as project_nonnegative, but this time the solution will be assigned to a new variable

function project_nonnegative(x::AbstractVector{T}) where T<:Real
    y = zeros(length(x))
    for i in eachindex(y)
        y[i] = max(zero(T), y[i])
    end
    return y
end

# Create the termination function
# which will take input A, b, c, x, y
# and will compute
# ϵ = ||Ax-b||/(1+||b||) + ||project_nonnegative(A'y-c)||/(1+||c||) + ||c'x - b'y||/(1+|c'x|+|b'y|)

function tolerance_LP(A::AbstractMatrix{T}, b::AbstractVector{T}, c::AbstractVector{T}, x::AbstractVector{T}, y::AbstractVector{T}) where T<:Real
    ϵ = norm(A*x-b,2)/(1+norm(b,2)) + norm(project_nonnegative(A'*y-c), 2)/(1+norm(c, 2)) + norm(c'*x - b'*y, 2)/(1+abs(c'*x)+abs(b'*y))
    return ϵ
end
