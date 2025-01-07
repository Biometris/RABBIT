using Revise
using MagicSimulate
using MagicBase
using Distributions
cd(@__DIR__)
pwd()

function get_breakpoints(homolog_fgl::AbstractVector)
    first.(homolog_fgl)
end

function fgl_seg_len(homolog_fgl::AbstractVector)
    diff([i[1] for i in homolog_fgl])
end


using Plots
using Plots.PlotMeasures

# matescheme = MateScheme(2,["Pairing","DH"],[1,1])
matescheme = MateScheme(2,["Pairing","Selfing"],[1,1])
matescheme = MateScheme(4,["Pairing","DH"],[2,1])
matescheme = MateScheme(8,["Pairing","Selfing"],[3,4])
# matescheme = MateScheme(8,["Pairing"],[3])
# matescheme = MateScheme(8,["Pairing","Sibling"],[2,10])
# matescheme = MateScheme(19,["FullDiallel", "RM1_E", "Selfing"],[1,4,6])
designcode = "4ril-self1"

pedinfo = matescheme
# pedinfo = designcode
@time contfgl,magicped = magicsimulate(pedinfo;
    popsize= 1000,        
    chrlen = 100*ones(5),
    outstem=nothing
)
# rm.(filter(x->occursin("sim",x), readdir()))

@time mosaic = plotcontfgl(contfgl;offspring=1:2)
@time savefig("mosaic.pdf")
display(mosaic)

using Plots
using LinearAlgebra

chr = 5
a = contfgl.offspringgeno[1][:,1:10]
breaks = get_breakpoints.(a)
n1 = length.(breaks)    
n2 = [length(unique(reduce(vcat,i))) for i in eachcol(breaks)]


nrecom1,nrecom2 = sum([begin     
    breaks = get_breakpoints.(a)
    n1 = reduce(vcat,length.(breaks)) .- 2    
    n2 = [length(unique(reduce(vcat,i))) for i in eachcol(breaks)] .- 2
    [n1,n2]
end for a in contfgl.offspringgeno])
xls1 = sort(unique(nrecom1))
yls1 = [1.0*sum(nrecom1 .== i) for i in xls1]
yls1 ./= 2
xls2 = sort(unique(nrecom2))
yls2 = [sum(nrecom2 .== i) for i in xls2]
gnrecom = plot([xls1,xls2],[yls1,yls2],
    # line = :scatter,    
    xlabel = "#recombination breakpoints",
    ylabel = "#offspring",
    label=["per haploid genome" "per diploid genome"],
)

segls = [begin 
    breaks = get_breakpoints.(a)
    seglen1 = reduce(vcat,diff.(breaks))
    seglen2 = reduce(vcat,[diff(sort(unique(reduce(vcat,i)))) for i in eachcol(breaks)])
    [seglen1,seglen2]    
end for a in contfgl.offspringgeno]
segls1 = reduce(vcat,first.(segls))
segls2 = reduce(vcat,last.(segls))
plot([segls1,segls2], seriestype=:stephist)
gseg = histogram([segls1,segls2],
    xlabel = "segement length (cM)",
    ylabel = "frequency",
    label=["along a single homolog" "along a pair of homologs"],
)

plot(mosaic,plot(gnrecom,gseg),layout=(2,1))