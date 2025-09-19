using Distributed
nprocs() < 11 && addprocs(11-nprocs())
@info string("nworkers=", nworkers())
@everywhere  using MagicCall


using Revise
using Distributions
using MagicBase
using MagicCall
cd(@__DIR__)
pwd()

using MagicReconstruct


using Plots
using StatsBase

ferr = 0.005
fixepso = 0.005
likeparasls = [    
    LikeParam(foundererror=ferr, offspringerror=fixepso,allelicbias=0.5, allelicoverdispersion=0.0,baseerror=0.001),        
    LikeParam(foundererror=ferr, offspringerror=fixepso,allelicoverdispersion=0.0,baseerror=0.001),  
    LikeParam(foundererror=ferr, offspringerror=fixepso,baseerror=0.001),
    # LikeParam(foundererror=ferr, offspringerror=fixepso,baseerror=nothing),
    # LikeParam(foundererror=ferr, offspringerror=nothing,baseerror=nothing)
]

function plotacc(resaccls, accnames)
    yls = [acc[accnames[1]][Symbol(accnames[2])] for acc in resaccls]
    xls = collect(eachindex(yls))
    xtickls = string.("M",eachindex(xls) .- 2)
    xtickls[1] = "Raw"
    plot(xls, yls; 
        marker=(5,:circle),     
        xlabel="LikeModel", 
        ylabel= join(accnames,"/"), 
        xticks=(xls, xtickls),
        legend= false,                
        yrange = (minimum(yls) * 0.0, maximum(yls)*1.05),
    )
end

# isfounderinbred = false

dataid = "sim"
genofile=string(dataid,"_magicsimulate_geno.vcf.gz")
pedfile = string(dataid,"_magicsimulate_ped.csv")
outstem = dataid*"_output"

resaccls = []

for likeindex in 1:(1+length(likeparasls))
    outstem = string(dataid,"_output_M", likeindex-1)
    israwcall = likeindex == 1
    @time magiccall(genofile,pedfile;
        isfounderinbred,              
        israwcall,
        isinfererror = !israwcall, 
        likeparam = likeparasls[israwcall ? 1 : likeindex-1],        
        # samplesize = 500,
        # burnin = 100,
        outstem,
    );    
    truegeno = formmagicgeno(dataid*"_magicsimulate_truegeno.csv.gz",pedfile; isfounderinbred);
    magicgeno = formmagicgeno(outstem*"_magiccall_geno.vcf.gz",pedfile; isfounderinbred, formatpriority=["GT"]);
    acc = magicaccuracy!(truegeno, magicgeno;alignfounder=true,isfounderinbred)
    push!(resaccls, acc)
    println("LikeModel = M",likeindex-2, ", imputation accuracy = ", acc)
    figerr = ploterrorcompare(truegeno, magicgeno)
    isnothing(figerr) || savefig(outstem*"_fig_error.png")
end

# mapfile = string(outstem, "_M2_magiccall_geno.vcf.gz")
# plotmarkererror(mapfile)


using Plots
namels = ["founder"=>"genoerr", "founder"=>"genomiss", "offspring"=>"genoerr", "offspring"=>"genomiss", "marker"=>"nungroup"]
gls = [plotacc(resaccls, name) for name in namels]
fig=plot(gls...;
    layout = (3,2), 
    size = (1000,1000),         
    legend = false,
    left_margin = 10Plots.mm,
    bottom_margin = 10Plots.mm,
)

outid = string("res_errormodel_depth",depth, "_eps",epsf, 
    isfounderinbred ? "_inbred" : "_outbred", "_baseerror", baseerror)
savefig(outid*"_acc_out.png")
fig

######################################################

# cd(@__DIR__)
# dataid = "sim"
# rm.(filter(x->occursin(dataid, x), readdir()))
# rm.(filter(x->occursin("outstem", x), readdir()))

