module HeatTransfer

export BC, solve_heat, solve_heat!, L_times!

using LinearAlgebra: I, mul!
using LinearMaps
import IterativeSolvers


"""
    BC

A boundary condition specified by a fixed temperature as a function of time.
"""
struct BC
    func::Function
end

"""
    BC(func)

Specify the boundary condition with a temperature function of time `func(t)`.

    BC(u::Real)

Specify the boundary condition as a constant temperature `u`.
"""
BC(u::Real) = BC(_ -> u)


"""
    solve_heat(initial::AbstractVector; len, t, bc, a=1, every=1)::Matrix

Solve the heat equation for an `initial` temperature state. Returns a matrix with
temperature states at each timestep in each column.

# Arguments
- `initial::AbstractVector`: Vector of initial temperatures.
- `len`: The length of the 1-D material.
- `t::AbstractRange`: Range of times to solve for, including the initial state.
- `bc::NTuple{2, BC}`: Boundary conditions at the edges of the 1-D material.
- `a=1`: Thermal diffusivity of the material.
- `every=1`: Return one state for every `every` states in `t`. For example, `every=2` saves
    every other state that is calculated.
"""
function solve_heat(initial::AbstractVector{T}; len::Real, t::AbstractRange{T},
                    bc::NTuple{2, BC}, a::Real=1, every::Int=1) where T<:AbstractFloat
    states = Matrix{T}(undef, length(initial), cld(length(t), every))
    return solve_heat!(states, initial, len, t, bc, a, every)
end

"""
    solve_heat(initial, x::AbstractRange; t, bc, a=1, every=1)::Matrix

`initial` can be provided as a function of position `initial(s)` with `x` as a range of
positions to use for calculations.
"""
function solve_heat(initial, x::AbstractRange{T}; kwargs...) where T
    x_incr = step(x) > 0 ? x : reverse(x)
    len = last(x_incr) - first(x_incr)
    return solve_heat(map(initial, x_incr); len=len, kwargs...)
end

"""
    solve_heat!(states::AbstractMatrix, initial::AbstractVector, len, t, bc, a=1, every=1)

[`solve_heat`](@ref) but store results in-place to `states`.
"""
function solve_heat!(states::AbstractMatrix{T}, initial::AbstractVector{T}, len::Real,
                     t::AbstractRange{T}, bc::NTuple{2, BC}, a::Real=1, every::Int=1
                    ) where T
    if size(states, 1) != length(initial)
        throw(DimensionMismatch("states must have the same amount of rows as initial."))
    end
    if size(states, 2) == 0
        return states
    end

    # Set the first column as the initial state
    copyto!(states, firstindex(states), initial, firstindex(initial), length(initial))

    # An interator over boundary condition temperatures (u1, u2) for each time
    current_bc, bc_iter = Iterators.peel(Tuple(b.func(t_k) for b in bc) for t_k in t)

    columns = Iterators.Stateful(eachcol(states))

    col = popfirst!(columns)
    col[begin], col[end] = current_bc

    u = col[begin+1:end-1]  # Allocate and initialize the inner temperatures

    Δx = len / (length(initial) - 1)
    k = a * step(t) / (2 * Δx^2)

    n = length(u)
    L = LinearMap(L_times!, n)
    A = I - k * L
    B = I + k * L
    y = Vector{T}(undef, n)

    for (i, next_bc) in enumerate(bc_iter)
        mul!(y, B, u)
        y[begin] += k * (current_bc[1] + next_bc[1])
        y[end] += k * (current_bc[2] + next_bc[2])

        IterativeSolvers.cg!(u, A, y)  # TODO: Add options for tolerance, maxiter, etc

        if i % every == 0
            col = popfirst!(columns)
            col[begin], col[end] = next_bc
            copyto!(col, firstindex(col) + 1, u, firstindex(u), n) 
        end

        current_bc = next_bc
    end

    return states
end

function L_times!(out::AbstractVector{T}, x::AbstractVector{T}) where T
    @assert eachindex(out) == eachindex(x) && length(x) > 1

    out[begin] = -2x[begin] + x[begin+1]
    for i in firstindex(x)+1 : lastindex(x)-1
        out[i] = x[i-1] - 2x[i] + x[i+1]
    end
    out[end] = x[end-1] - 2x[end]

    return out
end

end # module
