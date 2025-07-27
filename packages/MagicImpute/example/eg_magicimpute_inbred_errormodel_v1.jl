using Distributed
nprocs() < 6 && addprocs(6-nprocs())
@info string("nprocs=", nprocs())
@everywhere  using MagicReconstruct, MagicImpute

using Revise
using MagicBase, MagicReconstruct, MagicImpute
using StatsBase
using Distributions
cd(@__DIR__)
pwd()

likeparasls = [    
    LikeParam(offspringerror=0.005,allelicbias=0.5, allelicoverdispersion=0.0,baseerror=0.001),    
    LikeParam(allelicbias=0.5, allelicoverdispersion=0.0,baseerror=0.001),        
    LikeParam(allelicoverdispersion=0.0,baseerror=0.001),  
    LikeParam(baseerror=0.001),
    LikeParam(),                          
    # LikeParam(peroffspringerror = nothing)
]

function plotacc(resaccls, accnames)
    yls = [acc[accnames[1]][Symbol(accnames[2])] for acc in resaccls]
    xls = collect(eachindex(yls))
    plot(xls, yls; 
        marker=(5,:circle),     
        xlabel="LikeModel", 
        ylabel= join(accnames,"/"), 
        xticks=(xls, string.("M",eachindex(xls) .- 1)),
        legend= false,        
        yrange = (minimum(yls) * 0.9,maximum(yls)*1.05+ 0.01),
    )
end

isfounderinbred = false

dataid = "sim"
genofile=string(dataid,"_magicsimulate_geno.vcf.gz")
pedfile = string(dataid,"_magicsimulate_ped.csv")
outstem = dataid*"_output"


resaccls = []
for likeindex in 1:5
    outstem = string(dataid,"_output_M", likeindex-1)
    @time magicimpute(genofile,pedfile;
        isfounderinbred,      
        likeparam = likeparasls[likeindex],
        outstem,
    );    
    truegeno = formmagicgeno(dataid*"_magicsimulate_truegeno.csv.gz",pedfile; isfounderinbred);
    magicgeno = formmagicgeno(outstem*"_magicimpute_geno.vcf.gz",pedfile; isfounderinbred,formatpriority=["GT"]);
    acc = magicaccuracy!(truegeno, magicgeno;alignfounder=true,isfounderinbred)
    push!(resaccls, acc)
    println("LikeModel = M",likeindex-1, ", imputation accuracy = ", acc)
    figerr = ploterrorcompare(truegeno, magicgeno)
    isnothing(figerr) || savefig(outstem*"_fig_error.png")
end

resaccls[1]

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


# clear up

# cd(@__DIR__)
# dataid = "sim"
# rm.(filter(x->occursin(dataid, x), readdir()))

# dataid = "outstem"
# rm.(filter(x->occursin(dataid, x), readdir()))

