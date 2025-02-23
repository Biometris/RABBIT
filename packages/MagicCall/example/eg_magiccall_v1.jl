
# using Distributed
# nprocs() < 7 && addprocs(7-nprocs())
# @info string("nworkers=", nworkers())
# @everywhere  using MagicCall


using Revise
using MagicBase
using MagicCall
cd(@__DIR__)
pwd()

isfounderinbred = true
dataid = "sim"
genofile=string(dataid,"_magicfilter_geno.vcf.gz")
pedinfo = string(dataid, "_ped.csv")

@time magiccall(genofile,pedinfo;     
    isfounderinbred,
    likeparameters = LikeParameters(0.005, 0.04, 0.0, nothing, nothing, nothing, 0.0),
    model = "depmodel", 
)


calledgenofile = "outstem_magiccall_geno.vcf.gz"
magicgeno = formmagicgeno(calledgenofile,pedinfo;
    isfounderinbred, 
    formatpriority = ["GT"],
);
truefile = "sim_truegeno.csv.gz"
truegeno = formmagicgeno(truefile,pedinfo;isfounderinbred);

acc = magicaccuracy!(truegeno, magicgeno; isfounderinbred)
println(acc)

using Plots
using StatsBase
gls = [begin 
    trueerrls = copy(truegeno.markermap[1][!,col])    
    esterrls = copy(magicgeno.markermap[1][!,col])
    # if col == :offspringerror
    #    b = esterrls .> 0.5
    #    deleteat!(trueerrls,b)
    #    deleteat!(esterrls,b)
    # end
    r = round(cor(trueerrls,esterrls),digits=3)
    mtrue = round(mean(trueerrls),digits=3)
    mest = round(mean(esterrls),digits=3)
    @info string("col=", col, ",r=", r, ",mtrue=",mtrue,",mest=",mest)
    g = scatter(trueerrls,esterrls;
        xlab = string("true ", col),
        ylab = string("est ", col),
        label=string("r=", r),
    )
    g
end for col in [:offspringerror,:seqerror,:allelebalancemean,:allelebalancedisperse,:alleledropout]]
fig = plot(gls...,layout=(length(gls),1),
    size=(350,700),
    leftmargin = 5Plots.mm
)
display(fig)

