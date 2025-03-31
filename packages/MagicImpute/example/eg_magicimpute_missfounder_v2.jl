using Distributed
nprocs() < 7 && addprocs(7-nprocs())
@info string("nprocs=", nprocs())
@everywhere  using MagicReconstruct, MagicImpute

using Revise
using MagicBase, MagicReconstruct, MagicImpute
cd(@__DIR__)
pwd()

isfounderinbred = true


dataid = "si2"
genofile=string(dataid,"_magicsimulate_geno.vcf.gz")
pedfile = string(dataid,"_magicsimulate_ped.csv")
outstem = dataid*"_output"


magicgeno =formmagicgeno(genofile,pedfile; isfounderinbred); 
missingcode = isfounderinbred ? "N" : "NN"
for chr in eachindex(magicgeno.markermap)
    # nmiss = min(13, size(magicgeno.magicped.founderinfo,1))
    magicgeno.foundergeno[chr] .= missingcode
end
@time magicgeno = magicimpute!(magicgeno; 
    isfounderinbred,               
    # model = "depmodel",
    # model = ["depmodel","jointmodel"], 
    # byfounder = 8, 
    isallowmissing = false, 
    isrepeatimpute = true,    
    nrepeatmin = 6,
    nrepeatmax = 6,    
    isparallel = true, 
    isparallelfounder = true,       
    outstem, 
);


# accuracy
magicgeno = formmagicgeno(outstem*"_geno.vcf.gz",pedfile; isfounderinbred);
truegeno = formmagicgeno(dataid*"_magicsimulate_truegeno.csv.gz",pedfile; isfounderinbred);
facc = imputeaccuracy!(truegeno, magicgeno;isfounderinbred, alignfounder=true,targetfounder=true)
offacc = imputeaccuracy!(truegeno, magicgeno;isfounderinbred, alignfounder=true,targetfounder=false)
acc = magicaccuracy!(truegeno, magicgeno;alignfounder=true,isfounderinbred)
println(acc)
show(facc)
show(offacc)

using Plots
using StatsBase
colls = [:offspringerror,:seqerror,:allelebalancemean,:allelebalancedisperse,:alleledropout]
colls = colls[1:1]
gls = [begin 
    trueerrls = copy(truegeno.markermap[1][!,col])    
    esterrls = copy(magicgeno.markermap[1][!,col])
    b = .!ismissing.(trueerrls)
    if col == :alleledropout
        # b .= b .&& esterrls .> 1e-3    
        println("col=",col, ",#marker_incl=",sum(b))
    end    
    if any(b)
        r = cor(trueerrls[b],esterrls[b])
        mtrue = round(mean(trueerrls[b]),digits=3)
        mest = round(mean(esterrls[b]),digits=3)
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
display(fig)
savefig(outstem*"_fig_error.png")


genofile2 = outstem*"_geno.vcf.gz"
@time magicancestry = magicreconstruct(genofile2,pedfile;     
    # isignorephase = true,         
    isfounderinbred,
    nplot_subpop = 1, 
    outstem
);

truefgl = formmagicgeno(string(dataid,"_magicsimulate_truefgl.csv.gz"),pedfile; isfounderinbred);

acc2 = magicaccuracy!(truefgl, magicancestry;alignfounder=true,isfounderinbred)
println(acc2)


# MagicReconstruct.save_posterior_recom(magicancestry; isperchrom=true)

# using Plots
# noff = size(magicancestry.magicped.offspringinfo,1)
# offspring = rand(noff)
# plotcondprob(magicancestry; truefgl, probtype="diploprob",offspring=3)
# plot!(size=(1200,600))


# clear up

# cd(@__DIR__)
# dataid = "sim"
# rm.(filter(x->occursin(dataid, x), readdir()))

