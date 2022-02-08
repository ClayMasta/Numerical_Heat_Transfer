```@docs
BC(u::Real)
solve_heat(initial::AbstractVector{T}; len::Real, t::AbstractRange{T}, bc::NTuple{2, BC}, a::Real=1, every::Int=1) where T<:AbstractFloat
solve_heat(initial, x::AbstractRange{T}; kwargs...) where T
solve_heat!(states::AbstractMatrix{T}, initial::AbstractVector{T}, len::Real, t::AbstractRange{T}, bc::NTuple{2, BC}, a::Real=1, every::Int=1) where T
L_times!(out::AbstractVector{T}, x::AbstractVector{T}) where T
```