using Revise
using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
using SpectralOrdering
cd(@__DIR__)
pwd()
using Distributions, LinearAlgebra, StatsBase
using Plots

# Ctrl-S pauses the REPL terminal, so you can either put that keybinding
# into the whitelisted keybindings setting (under julia-client) or
# press Ctrl-Q to resume.

x = rand(Normal(0,1),1000)
# x=[1,2,3]
similarity = [exp(-(i-j)^2) for i=x, j=x]

trueorder= sortperm(x)
@time order1 = spectralordering(similarity;laplaciantype = "unnormalized")
@time order2 = spectralordering(similarity;laplaciantype = "randomwalk")
@time order3 = spectralordering(similarity;laplaciantype = "symmetric")

@time g1 = heatmap(similarity[order1, order1])
@time g2 = heatmap(similarity[order2, order2])
@time g3 = heatmap(similarity[order3, order3])
@time gt = heatmap(similarity[trueorder, trueorder])
@time plot(g1, g2, g3, gt)

plot(plot(x[order1],title="unnormalized"),
    plot(x[order2],title="randomwalk"),
    plot(x[order3],title="symmetric"),
    plot(x[trueorder],title="true"))

# println([corkendall(trueorder,j) for j=[order1,order2,order3]])
# println([corkendall(trueorder,reverse(j)) for j=[order1,order2,order3]])
