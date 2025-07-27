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

using Plots
using StatsBase

fixepso = 0.005

likeparasls = [    
    LikeParam(offspringerror=fixepso,allelicbias=0.5, allelicoverdispersion=0.0,baseerror=0.001),        
    LikeParam(offspringerror=fixepso,allelicoverdispersion=0.0,baseerror=0.001),  
    LikeParam(offspringerror=fixepso,baseerror=0.001),
    LikeParam(offspringerror=fixepso),
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
        yrange = (minimum(yls) * 0.5, maximum(yls)*1.05+ 0.01),
    )
end

isfounderinbred = false

dataid = "sim"
genofile=string(dataid,"_magicsimulate_geno.vcf.gz")
pedfile = string(dataid,"_magicsimulate_ped.csv")
outstem = dataid*"_output"

resaccls = []
threshoffspring = 0.9
for likeindex in eachindex(likeparasls)
    outstem = string(dataid,"_output_M", likeindex-1)
    @time magiccall(genofile,pedfile;
        isfounderinbred,              
        likeparam = likeparasls[likeindex],        
        threshoffspring, 
        outstem,
    );    
    truegeno = formmagicgeno(dataid*"_magicsimulate_truegeno.csv.gz",pedfile; isfounderinbred);
    magicgeno = formmagicgeno(outstem*"_magiccall_geno.vcf.gz",pedfile; isfounderinbred, formatpriority=["GT"]);
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
    isfounderinbred ? "_inbred" : "_outbred", "_fixepso",fixepso, "_thresh",threshoffspring, "_baseerror", baseerror)
savefig(outid*"_acc_out.png")
fig

######################################################


# cd(@__DIR__)
# dataid = "sim"
# rm.(filter(x->occursin(dataid, x), readdir()))
# rm.(filter(x->occursin("outstem", x), readdir()))

