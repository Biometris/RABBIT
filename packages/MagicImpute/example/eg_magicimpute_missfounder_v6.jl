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


##############################################################

using MagicCall
# genotype calling
magiccall(genofile,pedfile;             
    isfounderinbred, 
    likeparam = LikeParam(baseerror=0.001), 
    israwcall=true,
    isdelmonomorphic = true, 
    minmaf = 0.05, 
    maxomiss = 0.9, 
    maxhetero = 0.75, 
    # threshcall = 0.75,
    outstem, 
    # workdir,
)


genofile2 = string(outstem, "_magiccall_geno.vcf.gz")
magicgeno =formmagicgeno(genofile2,pedfile; isfounderinbred); 
for chr in eachindex(magicgeno.markermap)
    # magicgeno.foundergeno[chr] .= missingcode
    # magicgeno.foundergeno[chr][:,1:1] .= missingcode
    ls = magicgeno.foundergeno[chr]
    formatvec = magicgeno.markermap[chr][!,:founderformat]
    for i in 1:size(ls,1)
        if formatvec[i] == "GT_unphased"
            ls[i,1] = "NN"
        else
            ls[i,1] = Int16[0,0]
        end
        # ls[i,2] = [0,0]
    end
end
genofile3  = "temp.vcf.gz"
savegenodata(genofile3,magicgeno; keepcomment=false, infocomment = false)    


##############################################################

magicgeno =formmagicgeno(genofile3,pedfile; isfounderinbred); 
magicgeno.foundergeno[1][1:10,:]
magicgeno.offspringgeno[1][1:10,:]

@time magicimpute!(magicgeno;
    isfounderinbred,         
    isrepeatimpute = false, 
    isallowmissing = false,
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

