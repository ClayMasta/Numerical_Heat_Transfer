using Documenter

include("../src/NumericalHeatTransfer.jl")

using .HeatTransfer

makedocs(sitename = "Numerical Heat Transfer Documentation")

deploydocs(
    repo = "github.com/ClayMasta/Numerical_Heat_Transfer.git"
)