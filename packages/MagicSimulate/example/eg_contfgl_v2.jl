using Revise
using MagicSimulate
using MagicBase
using Distributions
cd(@__DIR__)
pwd()

using Plots

# matescheme = MateScheme(2,["Pairing","DH"],[1,1])
matescheme = MateScheme(2,["Pairing","Selfing"],[1,1])
matescheme = MateScheme(4,["Pairing","DH"],[2,1])
matescheme = MateScheme(4,["Pairing","Selfing"],[2,4])
matescheme = MateScheme(8,["Pairing","Selfing"],[3,1])
# matescheme = MateScheme(8,["Pairing"],[3])
# matescheme = MateScheme(8,["Pairing","Sibling"],[2,10])
# matescheme = MateScheme(19,["FullDiallel", "RM1_E", "Selfing"],[1,4,6])
# designcode = "4ril-self1"

pedinfo = matescheme
# pedinfo = designcode
@time contfgl,magicped = magicsimulate(pedinfo;
    popsize= 200,        
    chrlen = 100*ones(5),
    outstem=nothing
)
# rm.(filter(x->occursin("sim",x), readdir()))
contfgl.offspringgeno[1][1,1]


chr=1
a = contfgl.offspringgeno[chr]

breaks = MagicBase.get_breakpoints.(a)
nrecom = [length(unique(reduce(vcat,i))) - 2 for i in eachcol(breaks)]
histogram(nrecom)

nrecom = reduce(vcat,length.(breaks) .- 2)
histogram(nrecom)

g,_ = plotmosaic(contfgl;offspring=1:5)

g ,_ = hist_mosaic(contfgl)
