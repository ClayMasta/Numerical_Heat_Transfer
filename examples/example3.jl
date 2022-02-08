include("../src/heat-transfer.jl")

using .HeatTransfer
using Plots
using Printf            #Used for printing constant sigfigs for timer

export @userplot, @recipe

@userplot heatplot
@recipe function f(hp::heatplot)
    x, y, i, t, a = hp.args
    title --> "Heat Equation"
    xaxis --> ("x", (0,1))
    yaxis --> ("Temperature", (-1, 1))
    seriestype --> :line
    linewidth --> 2
    legend --> :none
    dpi --> 300
    c --> cgrad(:thermal)
    line_z --> y
    annotations --> [(0.1, -0.70, "α = " * string(a)), 
                    (0.137, -0.85, "T = " * @sprintf("%1.3f", t[i]) * " s")]
    annotationfontsize --> 12
    x, y
end

data = solve_heat(x -> x*(1-x), domain; t=time, bc=(BC(t -> cos(t)), BC(t -> sin(t))), a=a)

anim = @animate for i = 1 : length(time)
    heatplot(domain, data[:,i], i, time, a)
end

gif(anim, "example3.gif", fps=30)