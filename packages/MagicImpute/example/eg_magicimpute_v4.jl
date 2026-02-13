
# using Distributed
# nprocs() < 11 && addprocs(11-nprocs())
# @info string("nprocs=", nprocs())
# @everywhere  using MagicReconstruct, MagicImpute

using MagicBase, MagicReconstruct, MagicImpute
cd(@__DIR__)
pwd()

isfounderinbred = false

dataid = "sim"
genofile=string(dataid,"_magicsimulate_geno.vcf.gz")
pedfile = string(dataid,"_magicsimulate_ped.csv")
outstem = dataid*"_output"

magicgeno =formmagicgeno(genofile,pedfile; isfounderinbred); 
# missingcode = isfounderinbred ? "N" : "NN"
# for chr in eachindex(magicgeno.markermap)   
#     for j in [1,3]
#         magicgeno.foundergeno[chr][:,j] .= [[0,0] for _ in 1:size(magicgeno.foundergeno[chr],1)]
#     end
# end
@time magicimpute!(magicgeno; 
    model = "depmodel", 
    isfounderinbred,     
    isspacemarker = false,     
    isallowmissing = false, 
    # isdelmono = false, 
    isdelmarker = false, 
    tukeyfence = Inf, 
    threshlikeparam = ThreshLikeParam(offspringerror=1.0, allelicbias=1.0, allelicoverdispersion=Inf),             
    # target = "founder", 
    isrepeatimpute = false, 
    outstem,
);


magicgeno.markermap[1]

# accuracy
# magicgeno = formmagicgeno(outstem*"_geno.vcf.gz",pedfile; isfounderinbred);
truegeno = formmagicgeno(dataid*"_magicsimulate_truegeno.csv.gz",pedfile; isfounderinbred);
facc = imputeaccuracy!(truegeno, magicgeno;isfounderinbred, alignfounder=true,targetfounder=true)
offacc = imputeaccuracy!(truegeno, magicgeno;isfounderinbred, alignfounder=true,targetfounder=false)
acc = magicaccuracy!(truegeno, magicgeno;alignfounder=true,isfounderinbred)
println(acc)
show(facc)
show(offacc)

# if 1 == 0
#     magicgeno =formmagicgeno(genofile,pedfile; isfounderinbred); 
#     missingcode = isfounderinbred ? "N" : "NN"
#     for chr in eachindex(magicgeno.markermap)   
#         magicgeno.foundergeno[chr][:,1:2] .= missingcode
#     end
# end
@time magicancestry=magicreconstruct!(magicgeno;
    isfounderinbred,            
    outstem,
);

# magicancestry = readmagicancestry("outstem_magicreconstruct_ancestry.csv.gz");
truefgl = formmagicgeno(string(dataid,"_magicsimulate_truefgl.csv.gz"),pedfile; isfounderinbred);
acc = magicaccuracy!(truefgl, magicancestry; isfounderinbred)
@info "" acc
plotcondprob(magicancestry;
    truefgl, 
    offspring = 1,
    probtype="haploprob"
)

# clear up
# cd(@__DIR__)
# dataid = "sim"
# rm.(filter(x->occursin(dataid, x), readdir()))
# rm.(filter(x->occursin("outstem", x), readdir()))


using CSV, DataFrames
using Plots
using MagicBase
using StatsBase
cd(@__DIR__)
pwd()
exportall(MagicBase)

popsize = 200
isfounderinbred = false
datadir = joinpath(@__DIR__, "simdata")
resdir = joinpath(@__DIR__, "res_magicimpute")
# resdir = joinpath(@__DIR__, "res_magicimpute_mma")



withod = 0
outstemls = [string("simgeno_isgbs1_withod",i,"_popsize", popsize,) for i in [0,1]]
pushfirst!(outstemls, string("simgeno_isgbs0_withod0_popsize", popsize))



alignfounder = true
outstem = outstemls[1]
truegenofile = joinpath(datadir, string(outstem, "_magicsimulate_truegeno.csv.gz"))   
genofile = joinpath(resdir, string(outstem, "_magicimpute_geno.vcf.gz"))
# genofile = joinpath(resdir, string(outstem, "_out_ImputedGenotype.vcf.gz"))
pedfile = joinpath(datadir, string(outstem, "_magicsimulate_ped.csv"))                        

magicgeno = formmagicgeno(genofile, pedfile; isfounderinbred,formatpriority=["GT"]);
MagicBase.rawgenocall!(magicgeno; targets = ["offspring"], callthreshold = 0.0, isfounderinbred)
truegeno = formmagicgeno(truegenofile, pedfile; isfounderinbred);
alignoffspring!(truegeno,magicgeno) 
alignchromosome!(truegeno,magicgeno)     
alignmarker!(truegeno,magicgeno)     
alignfounder && alignfounder!(truegeno,magicgeno)
merge_chromosome!(magicgeno)
merge_chromosome!(truegeno)

isfounder = false
chr = 1
if isfounder 
    chrest = join.(sort.(magicgeno.foundergeno[chr])) 
    chrtrue = join.(sort.(truegeno.foundergeno[chr]))
else
    chrest = magicgeno.offspringgeno[chr]
    chrtrue = join.(sort.(truegeno.offspringgeno[chr]))
end
gset = unique(chrest)
issubset(gset,["11","12","22","1N","2N","NN"]) || error("unexpected imputed genotypes: ",gset)

estls = eachrow(chrest)
truels = eachrow(chrtrue)  
accls = [begin 
    ls = map((x,y)->occursin("N",x) ? missing : x==y, estls[i], truels[i])
    mean(skipmissing(ls))
end for i in eachindex(estls,truels)]
ar2ls = [begin     
    infols = split(i,";")
    filter!(x->occursin("AR2=",x), infols)
    parse(Float64,replace(infols[1],"AR2="=>"") )
end for i in magicgeno.markermap[chr][!,:info]]
# ar2ls ./= missls    
b = isnan.(ar2ls)
if any(b)
    deleteat!(ar2ls,b) 
    deleteat!(accls,b) 
end
plot(ar2ls, accls, seriestype=:scatter, 
    # xlim=(minimum(ar2ls)*0.999,1.001), ylim=(minimum(accls)*0.999,1.001), 
    # yrange = (0.8,1.001),
    # yrange = (0.987,0.992),
    # yrange = (0.999,1.001),
    xlabel="AR2", ylabel="Imputation accuracy", 
    # xscale = :log10,
    # yscale = :log10,
) 
cor(ar2ls, accls)



