using Distributed
nprocs() < 6 && addprocs(6-nprocs())
@info string("nprocs=", nprocs())
@everywhere  using MagicReconstruct, MagicImpute

using Revise
using MagicBase, MagicReconstruct, MagicImpute
cd(@__DIR__)
pwd()

isfounderinbred = true


dataid = "sim"
genofile=string(dataid,"_magicsimulate_geno.vcf.gz")
pedfile = string(dataid,"_magicsimulate_ped.csv")
outstem = dataid*"_output"


magicgeno =formmagicgeno(genofile,pedfile; isfounderinbred); 
missingcode = isfounderinbred ? "N" : "NN"
for chr in eachindex(magicgeno.markermap)
    # magicgeno.foundergeno[chr] .= missingcode
    magicgeno.foundergeno[chr][:,1:9] .= missingcode
end
@time magicgeno = magicimpute!(magicgeno; 
    isfounderinbred,     
    # target = "founder",      
    # likeparameters = LikeParam(peroffspringerror=0),      
    model = "depmodel",    
    # isgreedy = true, 
    # byfounder = 4, 
    # threshproposal = 0.9, 
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

println("")
for i in eachindex(magicgeno.markermap)
    acc = magicaccuracy!(truegeno, magicgeno; alignfounder=true,chrsubset=[i])    
    println("chr=", i, ";", acc["founder"])
end
0

# clear up

# cd(@__DIR__)
# dataid = "sim"
# rm.(filter(x->occursin(dataid, x), readdir()))

