# using Distributed
# nprocs() < 3 && addprocs(3-nprocs())
# @info string("nprocs=", nprocs())
# @everywhere  using MagicReconstruct, MagicImpute

using Revise
using MagicBase, MagicReconstruct, MagicImpute
using StatsBase
cd(@__DIR__)
pwd()

isfounderinbred = false


dataid = "sim"
genofile=string(dataid,"_magicsimulate_geno.vcf.gz")
pedfile = string(dataid,"_magicsimulate_ped.csv")
outstem = dataid*"_output2"

tused = @elapsed magicimpute(genofile,pedfile;
    isfounderinbred,  
    # model = "depmodel", 
    # likeparameters = LikeParameters(peroffspringerror=0.0),              
    # isbinning = true, 
    # bincm = 0.05, 
    # isinfererror = false, 
    # byfounder = 2, 
    # isrepeatimpute = false, 
    # nrepeatmin = 3, 
    nrepeatmax = 1,     
    # iscorrectfounder = true, 
    # isdelmarker = false,     
    # tempdirectory = "D://Temp",    
    outstem,
)



# inputmagicgeno = formmagicgeno(genofile,pedfile;isfounderinbred); 
# magicgeno = deepcopy(inputmagicgeno);
# nparentmiss = 1
# nchr = length(magicgeno.markermap)
# if nparentmiss>0
#     nparent = size(magicgeno.magicped.founderinfo,1)
#     for chr = 1:nchr
#         jj = sample(1:nparent,nparentmiss,replace=false)
#         println("jj=",jj)
#         formatls = magicgeno.markermap[chr][!,:founderformat]
#         misscodels = MagicBase.get_missingcode(formatls)
#         magicgeno.foundergeno[chr][:,jj] .= misscodels                    
#     end
# end

# tused = @elapsed magicimpute!(magicgeno;    
#     isfounderinbred,            
#     # byfounder = 1, 
#     isrepeatimpute = false, 
#     # nrepeatmax = 4, 
#     # iscorrectfounder = false,
#     tempdirectory = "D://Temp",        
#     outstem, 
# )


# accuracy
magicgeno = formmagicgeno(outstem*"_magicimpute_geno.vcf.gz",pedfile; isfounderinbred);
truegeno = formmagicgeno(dataid*"_magicsimulate_truegeno.csv.gz",pedfile; isfounderinbred);
facc = imputeaccuracy!(truegeno, magicgeno;isfounderinbred, alignfounder=true,targetfounder=true)
offacc = imputeaccuracy!(truegeno, magicgeno;isfounderinbred, alignfounder=true,targetfounder=false)
acc = magicaccuracy!(truegeno, magicgeno;alignfounder=true,isfounderinbred)
println(acc)


using Plots
using StatsBase
colls = [:offspringerror,:seqerror,:allelebalancemean,:allelebalancedisperse,:alleledropout]
# colls = colls[1:1]
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



genofile2 = outstem*"_magicimpute_geno.vcf.gz"
@time magicancestry = magicreconstruct(genofile2,pedfile;     
    isignorephase = true,         
    isfounderinbred,
    nplot_subpop = 1, 
    outstem
);


using Plots
truefgl = formmagicgeno(string(dataid,"_magicsimulate_truefgl.csv.gz"),pedfile; isfounderinbred);
acc2 = magicaccuracy!(truefgl, magicancestry;alignfounder=true,isfounderinbred)
println(acc2)

noff = size(magicancestry.magicped.offspringinfo,1)
offspring = rand(1:noff)
offspring = 40-4
plotcondprob(magicancestry; truefgl, probtype="genoprob",offspring)
plot!(size=(1200,600))

# clear up

cd(@__DIR__)
dataid = "sim"
rm.(filter(x->occursin(dataid, x), readdir()))

