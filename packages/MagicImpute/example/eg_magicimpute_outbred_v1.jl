# using Distributed
# nprocs() < 3 && addprocs(3-nprocs())
# @info string("nprocs=", nprocs())
# @everywhere  using MagicReconstruct, MagicImpute

using Revise
using MagicBase, MagicReconstruct, MagicImpute
using StatsBase
using Distributions
cd(@__DIR__)
pwd()


like = "LikeParam(offspringerror=0.1)"
MagicBase.parse_likeparameters(like, "LikeParam")




isfounderinbred = false

dataid = "sim"
genofile=string(dataid,"_magicsimulate_geno.vcf.gz")
pedfile = string(dataid,"_magicsimulate_ped.csv")
outstem = dataid*"_output"

@time magicimpute(genofile,pedfile;
    isfounderinbred,      
    # model = "depmodel", 
    # likeparameters = LikeParam(foundererror=0.01),
    # likeparameters = LikeParam(peroffspringerror=0.0,allelicbias=0.5, allelicoverdispersion=0),            
    # priorlikeparameters = PriorLikeParam(offspringerror=Beta(1,99), baseerror=Beta(1,999)), 
    # threshlikeparameters = ThreshLikeParam(allelicoverdispersion=2),
    # isbinning = true, 
    # bincm = 0.05, 
    # isinfererror = false, 
    # byfounder = 2, 
    # isrepeatimpute = true, 
    # nrepeatmin = 2, 
    # nrepeatmax = 2,         
    # isspacemarker = true,
    # isordermarker = true,
    # iscorrectfounder = true, 
    # isdelmarker = false,     
    # tempdirectory = "D://Temp",    
    # isgreedy = true,     
    outstem,
);

# inputmagicgeno = formmagicgeno(genofile,pedfile;isfounderinbred); 
# magicgeno = deepcopy(inputmagicgeno);
# nparentmiss = 2
# nchr = length(magicgeno.markermap)
# if nparentmiss>0
#     nparent = size(magicgeno.magicped.founderinfo,1)
#     for chr = 1:nchr
#         jj = sample(1:nparent,nparentmiss,replace=false)
#         # println("jj=",jj)
#         formatls = magicgeno.markermap[chr][!,:founderformat]
#         misscodels = MagicBase.get_missingcode(formatls)
#         magicgeno.foundergeno[chr][:,jj] .= misscodels                    
#     end
# end

# tused = @elapsed magicimpute!(magicgeno;    
#     isfounderinbred,            
#     # byfounder = 1, 
#     isrepeatimpute = true, 
#     # nrepeatmax = 4, 
#     # iscorrectfounder = false,    
#     isordermarker = true,
#     isspacemarker = true,
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

savefig(outstem*"_fig_error.png")
display(fig)


# genofile2 = outstem*"_magicimpute_geno.vcf.gz"
# @time magicancestry = magicreconstruct(genofile2,pedfile;     
#     isignorephase = true,         
#     isfounderinbred,
#     nplot_subpop = 1, 
#     outstem
# );


# using Plots
# truefgl = formmagicgeno(string(dataid,"_magicsimulate_truefgl.csv.gz"),pedfile; isfounderinbred);
# acc2 = magicaccuracy!(truefgl, magicancestry;alignfounder=true,isfounderinbred)
# println(acc2)

# noff = size(magicancestry.magicped.offspringinfo,1)
# offspring = rand(1:noff)
# plotcondprob(magicancestry; truefgl, probtype="genoprob",offspring)
# plot!(size=(1200,600))



# clear up

# cd(@__DIR__)
# dataid = "sim"
# rm.(filter(x->occursin(dataid, x), readdir()))

# dataid = "outstem"
# rm.(filter(x->occursin(dataid, x), readdir()))

