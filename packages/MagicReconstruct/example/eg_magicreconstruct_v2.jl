# using Distributed
# nprocs() < 6 && addprocs(6-nprocs())
# @info string("nworkers=", nworkers())
# @everywhere  using MagicReconstruct

using Revise
using MagicBase, MagicReconstruct
cd(@__DIR__)
pwd()

dataid = "sim"
genofile=dataid*"_magicsimulate_geno.vcf.gz"
pedinfo = dataid*"_magicsimulate_ped.csv"
@time magicancestry=magicreconstruct(genofile,pedinfo;
    likeparameters = LikeParam(offspringerror=0.995),
    model = "depmodel",
    # isfounderinbred,        
    # hmmalg = "viterbi",    
);

# magicancestry = readmagicancestry("outstem_magicreconstruct_ancestry.csv.gz");
magicancestry = readmagicancestry("outstem_magicreconstruct_ancestry.csv.gz");
truefgl = formmagicgeno(string(dataid,"_magicsimulate_truefgl.csv.gz"),pedinfo);
acc = magicaccuracy!(truefgl, magicancestry)
println(acc)

plotcondprob(magicancestry;
    truefgl, 
    offspring=1,
    probtype="haploprob"
)


# cd(@__DIR__)
# dataid = "sim"
# rm.(filter(x->occursin(dataid,x), readdir()))
# rm.(filter(x->occursin("outstem",x), readdir()))
