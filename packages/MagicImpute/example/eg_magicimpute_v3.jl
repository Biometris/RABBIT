
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
# magicgeno.foundergeno[2]
# missingcode = isfounderinbred ? "N" : "NN"
# for chr in eachindex(magicgeno.markermap)   
#     for j in [1,3]
#         magicgeno.foundergeno[chr][:,j] .= [[0,0] for _ in 1:size(magicgeno.foundergeno[chr],1)]
#     end
# end
@time magicimpute!(magicgeno; 
    isfounderinbred,     
    isspacemarker = false,     
    isallowmissing = false, 
    isdelmono = false, 
    isdelmarker = false, 
    tukeyfence = Inf, 
    threshlikeparam = ThreshLikeParam(offspringerror=1.0, allelicbias=1.0, allelicoverdispersion=Inf),             
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

