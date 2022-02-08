include("../src/NumericalHeatTransfer.jl")

using .HeatTransfer
using Plots
using Printf

@userplot heatplot
@recipe function f!(hp::heatplot)
    x, y, i, t, a = hp.args
    title --> "Heat Equation"
    xaxis --> ("x", (-1,1))
    yaxis --> ("Temperature", (-0.1, 1))
    seriestype --> :line
    linewidth --> 2
    legend --> :none
    dpi --> 300
    c --> cgrad(:thermal)
    line_z --> y
    annotations --> [(-0.75, 0.1, "α = " * string(a)), 
                    (-0.725, 0, "T = " * @sprintf("%1.3f", t[i]) * " s")]
    annotationfontsize --> 12
    x, y
end

mean = 0.5
std = 0.25
domain = LinRange(-1,1,400)
time = range(0,7;step=0.025)
a = 0.125

data = solve_heat(x -> ( 1 / (std * sqrt(2 * π)) ) * exp(-0.5 * ((x-mean)/std)^2), 
                    domain; t=time, bc=(BC(0), BC(0)), a=a)

anim = @animate for i = 1 : length(time)
    heatplot(domain, data[:,i], i, time, a)
end

gif(anim, "example4.gif", fps=30)