

using Revise
using MagicBase, MagicReconstruct, MagicImpute
cd(@__DIR__)
pwd()

# isfounderinbred = true


dataid = "sim"
genofile=string(dataid,"_magicsimulate_geno.vcf.gz")
pedfile = string(dataid,"_magicsimulate_ped.csv")
outstem = dataid*"_output"

@time magicgeno = magicimpute(genofile,pedfile;    
    isfounderinbred,            
    target = "founder", 
    # model = "depmodel",
    # likeparameters = LikeParam(peroffspringerror=nothing),
    # priorlikeparameters = PriorLikeParam(offspringerror=Beta(1,19), 
    #     baseerror=Beta(1,999),allelicoverdispersion=Exponential(0.3)),
    # priorlikeparameters = PriorLikeParam(allelicoverdispersion=Exponential(0.5)),  
    isspacemarker = true, 
    spacebyviterbi = true, 
    # skeletonsize = 10^6,
    outstem,
);



# accuracy
# magicgeno = formmagicgeno(outstem*"_magicimpute_geno.vcf.gz",pedfile; isfounderinbred);
truegeno = formmagicgeno(dataid*"_magicsimulate_truegeno.csv.gz",pedfile; isfounderinbred);
facc = imputeaccuracy!(truegeno, magicgeno;isfounderinbred, alignfounder=true,targetfounder=true)
offacc = imputeaccuracy!(truegeno, magicgeno;isfounderinbred, alignfounder=true,targetfounder=false)
acc = magicaccuracy!(truegeno, magicgeno;alignfounder=true,isfounderinbred)
println(acc)
show(facc)
show(offacc)

# using Plots
# using StatsBase
# colls = [:offspringerror,:baseerror,:allelicbias,:allelicoverdispersion,:allelicdropout]
# # colls = colls[1:1]
# # peroffls = magicgeno.magicped.offspringinfo[!,:peroffspringerror_1] # only one LG
# gls = [begin 
#     trueerrls = copy(truegeno.markermap[1][!,col])    
#     esterrls = copy(magicgeno.markermap[1][!,col])
#     # if col == :offspringerror
#     #     esterrls .= [ismissing(i) ? i : mean(1 .- (1 .- peroffls) .* (1 - i)) for i in esterrls]
#     # end
#     b = .!ismissing.(trueerrls)
#     if col == :allelicdropout
#         # b .= b .&& esterrls .> 1e-3    
#         println("col=",col, ",#marker_incl=",sum(b))
#     end    
#     if any(b)
#         r = cor(trueerrls[b],esterrls[b])
#         mtrue = round(mean(trueerrls[b]),digits=4)
#         mest = round(mean(esterrls[b]),digits=4)
#         @info string("col=", col, ",mtrue=",mtrue,",mest=",mest)
#         y = esterrls[b] 
#         # y = esterrls[b] - trueerrls[b]
#         g = scatter(trueerrls[b],y;
#             xlab = string("true ", col),
#             ylab = string("est ", col),
#             label=string("r=", round(r,digits=3)),
#         )
#         g
#     else
#         nothing
#     end
# end for col in colls]
# deleteat!(gls, isnothing.(gls))
# fig = plot(gls...,layout=(length(gls),1),
#     size=(350,700),
#     leftmargin = 5Plots.mm
# )
# display(fig)
# savefig(outstem*"_fig_error.png")

# ls = magicgeno.magicped.offspringinfo[!,:peroffspringerror_1]
# histogram(ls)
# savefig(outstem*"_fig_error_peroff.png")

# clear up


cd(@__DIR__)
dataid = "sim"
rm.(filter(x->occursin(dataid, x), readdir()))

dataid = "outstem"
rm.(filter(x->occursin(dataid, x), readdir()))

