using Documenter

include("../src/heat-transfer.jl")

using .HeatTransfer

makedocs(sitename = "Numerical Heat Transfer Documentation")

deploydocs(
    repo = "github.com/arturofburgos/Julia_Repo_Studies.git"
)