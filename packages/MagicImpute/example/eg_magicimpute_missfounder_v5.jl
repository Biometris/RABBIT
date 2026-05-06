# using Distributed
# nprocs() < 6 && addprocs(6-nprocs())
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
outstem = dataid*"_output"


magicgeno =formmagicgeno(genofile,pedfile; isfounderinbred); 
missingcode = isfounderinbred ? "N" : "NN"
magicgeno.foundergeno[1][1,1:2]
for chr in eachindex(magicgeno.markermap)
    # magicgeno.foundergeno[chr] .= missingcode
    # magicgeno.foundergeno[chr][:,1:1] .= missingcode
    ls = magicgeno.foundergeno[chr]
    for i in 1:size(ls,1)
        ls[i,1] = [0,0]
        # ls[i,2] = [0,0]
    end
end

using MagicCall
genofile2  = "temp.vcf.gz"
savegenodata(genofile2,magicgeno; keepcomment=false, infocomment = false)    
# genotype calling
magiccall(genofile2,pedfile;             
    isfounderinbred, 
    likeparam = LikeParam(baseerror=0.001), 
    israwcall=true,
    isdelmonomorphic = true, 
    minmaf = 0.05, 
    maxomiss = 0.95, 
    maxhetero = 0.7, 
    outstem, 
    # workdir,
)
genofile3 = string(outstem, "_magiccall_geno.vcf.gz")

magicgeno =formmagicgeno(genofile3,pedfile; isfounderinbred); 
magicgeno.foundergeno[1][1:10,:]
magicgeno.offspringgeno[1][1:10,:]

@time magicimpute!(magicgeno;
    isfounderinbred,     
    isrepeatimpute = false, 
    blockdelmarker = true,     
    # target = "founder",          
    # model = "depmodel",        
    outstem, 
);
0

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

