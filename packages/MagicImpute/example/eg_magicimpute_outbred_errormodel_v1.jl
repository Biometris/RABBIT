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

using Plots
using StatsBase

likeparasls = [    
    LikeParam(offspringerror=0.005,allelicbias=0.5, allelicoverdispersion=0.0,baseerror=0.001),    
    LikeParam(allelicbias=0.5, allelicoverdispersion=0.0,baseerror=0.001),        
    LikeParam(allelicoverdispersion=0.0,baseerror=0.001),  
    LikeParam(baseerror=0.001),
    LikeParam(),                          
    LikeParam(peroffspringerror = nothing)
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
    )
end

function ploterracc(magicgeno,truegeno)
    colls = [:offspringerror,:baseerror,:allelicbias,:allelicoverdispersion,:allelicdropout]
    # colls = colls[1:1]
    gls = [begin 
        trueerrls = copy(truegeno.markermap[1][!,col])    
        esterrls = copy(magicgeno.markermap[1][!,col])
        b = .!ismissing.(trueerrls)
        if col == :allelicdropout
            # b .= b .&& esterrls .> 1e-3    
            println("col=",col, ",#marker_incl=",sum(b))
        end    
        if any(b)
            r = cor(trueerrls[b],esterrls[b])
            mtrue = round(mean(trueerrls[b]),digits=4)
            mest = round(mean(esterrls[b]),digits=4)
            @info string("col=", col, ",mtrue=",mtrue,",mest=",mest)
            g = scatter(trueerrls[b],esterrls[b];
                xlab = string("true ", col),
                ylab = string("est ", col),
                label=string("r=", round(r,digits=3)),
            )
            g
        else
            nothing
        end
    end for col in colls]
    deleteat!(gls, isnothing.(gls))
    fig = plot(gls...,layout=(length(gls),1),
        size=(350,700),
        leftmargin = 5Plots.mm
    )
    fig
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
        # isrepeatimpute = true, 
        likeparam = likeparasls[likeindex],
        # softthreshlikeparameters = SoftThreshLikeParam(allelicbias=0.75),
        # threshlikeparameters = ThreshLikeParam(offspringerror=0.25,baseerror=0.25,allelicbias=0.95,allelicoverdispersion=2),
        outstem,
    );    
    truegeno = formmagicgeno(dataid*"_magicsimulate_truegeno.csv.gz",pedfile; isfounderinbred);
    magicgeno = formmagicgeno(outstem*"_magicimpute_geno.vcf.gz",pedfile; isfounderinbred,formatpriority=["GT"]);
    acc = magicaccuracy!(truegeno, magicgeno;alignfounder=true,isfounderinbred)
    push!(resaccls, acc)
    println("LikeModel = M",likeindex-1, ", imputation accuracy = ", acc)
    figerr = ploterracc(magicgeno,truegeno)
    savefig(outstem*"_fig_error.png")
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
# clear up

# cd(@__DIR__)
# dataid = "sim"
# rm.(filter(x->occursin(dataid, x), readdir()))

# dataid = "outstem"
# rm.(filter(x->occursin(dataid, x), readdir()))


