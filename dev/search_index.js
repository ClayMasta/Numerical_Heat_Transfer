var documenterSearchIndex = {"docs":
[{"location":"","page":"-","title":"-","text":"BC(u::Real)\nsolve_heat(initial::AbstractVector{T}; len::Real, t::AbstractRange{T}, bc::NTuple{2, BC}, a::Real=1, every::Int=1) where T<:AbstractFloat\nsolve_heat(initial, x::AbstractRange{T}; kwargs...) where T\nsolve_heat!(states::AbstractMatrix{T}, initial::AbstractVector{T}, len::Real, t::AbstractRange{T}, bc::NTuple{2, BC}, a::Real=1, every::Int=1) where T\nL_times!(out::AbstractVector{T}, x::AbstractVector{T}) where T","category":"page"},{"location":"#Main.HeatTransfer.BC-Tuple{Real}","page":"-","title":"Main.HeatTransfer.BC","text":"BC(func)\n\nSpecify the boundary condition with a temperature function of time func(t).\n\nBC(u::Real)\n\nSpecify the boundary condition as a constant temperature u.\n\n\n\n\n\n","category":"method"},{"location":"#Main.HeatTransfer.solve_heat-Union{Tuple{AbstractVector{T}}, Tuple{T}} where T<:AbstractFloat","page":"-","title":"Main.HeatTransfer.solve_heat","text":"solve_heat(initial::AbstractVector; len, t, bc, a=1, every=1)::Matrix\n\nSolve the heat equation for an initial temperature state. Returns a matrix with temperature states at each timestep in each column.\n\nArguments\n\ninitial::AbstractVector: Vector of initial temperatures.\nlen: The length of the 1-D material.\nt::AbstractRange: Range of times to solve for, including the initial state.\nbc::NTuple{2, BC}: Boundary conditions at the edges of the 1-D material.\na=1: Thermal diffusivity of the material.\nevery=1: Return one state for every every states in t. For example, every=2 saves   every other state that is calculated.\n\n\n\n\n\n","category":"method"},{"location":"#Main.HeatTransfer.solve_heat-Union{Tuple{T}, Tuple{Any, AbstractRange{T}}} where T","page":"-","title":"Main.HeatTransfer.solve_heat","text":"solve_heat(initial, x::AbstractRange; t, bc, a=1, every=1)::Matrix\n\ninitial can be provided as a function of position initial(s) with x as a range of positions to use for calculations.\n\n\n\n\n\n","category":"method"},{"location":"#Main.HeatTransfer.solve_heat!-Union{Tuple{T}, Tuple{AbstractMatrix{T}, AbstractVector{T}, Real, AbstractRange{T}, Tuple{BC, BC}}, Tuple{AbstractMatrix{T}, AbstractVector{T}, Real, AbstractRange{T}, Tuple{BC, BC}, Real}, Tuple{AbstractMatrix{T}, AbstractVector{T}, Real, AbstractRange{T}, Tuple{BC, BC}, Real, Int64}} where T","page":"-","title":"Main.HeatTransfer.solve_heat!","text":"solve_heat!(states::AbstractMatrix, initial::AbstractVector, len, t, bc, a=1, every=1)\n\nsolve_heat but store results in-place to states.\n\n\n\n\n\n","category":"method"}]
}
