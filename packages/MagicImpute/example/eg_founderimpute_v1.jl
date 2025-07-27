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
    # target = "founder",
    # model = "depmodel", 
    # likeparameters = LikeParam(peroffspringerror=0.0),              
    # isbinning = true, 
    # bincm = 0.05, 
    # isinfererror = false, 
    # byfounder = 2, 
    isrepeatimpute = false, 
    # nrepeatmin = 3, 
    nrepeatmax = 5,     
    # iscorrectfounder = true, 
    # isdelmarker = false,     
    # tempdirectory = "D://Temp",    
    outstem,
)



# accuracy
magicgeno = formmagicgeno(outstem*"_magicimpute_geno.vcf.gz",pedfile; isfounderinbred);
truegeno = formmagicgeno(dataid*"_magicsimulate_truegeno.csv.gz",pedfile; isfounderinbred);
facc = imputeaccuracy!(truegeno, magicgeno;isfounderinbred, alignfounder=true,targetfounder=true)
offacc = imputeaccuracy!(truegeno, magicgeno;isfounderinbred, alignfounder=true,targetfounder=false)
show(facc)
show(offacc)
acc = magicaccuracy!(truegeno, magicgeno;alignfounder=true,isfounderinbred)
show(acc)