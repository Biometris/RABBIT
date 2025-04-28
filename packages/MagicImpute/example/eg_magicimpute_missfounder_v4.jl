# using Distributed
# nprocs() < 7 && addprocs(7-nprocs())
# @info string("nprocs=", nprocs())
# @everywhere  using MagicReconstruct, MagicImpute

using Revise
using MagicBase, MagicReconstruct, MagicImpute
cd(@__DIR__)
pwd()

isfounderinbred = false


dataid = "sim"
genofile=string(dataid,"_magicsimulate_geno.vcf.gz")
pedfile = string(dataid,"_magicsimulate_ped.csv")
outstem = dataid*"_output2"


magicgeno =formmagicgeno(genofile,pedfile; isfounderinbred); 
missingcode = isfounderinbred ? "N" : "NN"
for chr in eachindex(magicgeno.markermap)
    # nmiss = min(13, size(magicgeno.magicped.founderinfo,1))
    magicgeno.foundergeno[chr].= missingcode
end
@time magicgeno = magicimpute!(magicgeno; 
    isfounderinbred,     
    # target = "founder",           
    # model = "depmodel",    
    # byfounder = 2, 
    # threshproposal = 0.55, 
    # isallowmissing = true,     
    # isrepeatimpute = true, 
    # startbyhalf = 2,     
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


# clear up

# cd(@__DIR__)
# dataid = "sim"
# rm.(filter(x->occursin(dataid, x), readdir()))

