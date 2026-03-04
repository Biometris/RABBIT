

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
ingeno = deepcopy(magicgeno);

missingcode = isfounderinbred ? "N" : "NN"
for chr in eachindex(magicgeno.markermap)   
    magicgeno.foundergeno[chr][:,1:1] .= missingcode
end
@time MagicImpute.magicimpute!(magicgeno; 
    isfounderinbred,           
    isspacemarker = false,     
    isallowmissing = false, 
    isrepeatimpute = false, 
    isdelmono = false, 
    isdelmarker = false, 
    tukeyfence = Inf, 
    threshlikeparam = ThreshLikeParam(offspringerror=1.0, allelicbias=1.0, allelicoverdispersion=Inf),    
    threshar2 = 0.0,        
    iscorrectfounder = true,   
    # threshproposal = 0.4,
    # byfounder = 2, 
    outstem,
);
0

# accuracy
# magicgeno = formmagicgeno(outstem*"_magicimpute_geno.vcf.gz",pedfile; isfounderinbred);
truegeno = formmagicgeno(dataid*"_magicsimulate_truegeno.csv.gz",pedfile; isfounderinbred);
facc = imputeaccuracy!(truegeno, magicgeno;isfounderinbred, alignfounder=true,targetfounder=true)
offacc = imputeaccuracy!(truegeno, magicgeno;isfounderinbred, alignfounder=true,targetfounder=false)
acc = magicaccuracy!(truegeno, magicgeno;alignfounder=true,isfounderinbred)
println(acc)
show(facc)
show(offacc)


genofile2 = outstem*"_geno.vcf.gz"
@time magicancestry=magicreconstruct(genofile,pedfile;
    isfounderinbred,
    # model = "depmodel",    
);

# magicancestry = readmagicancestry("outstem_magicreconstruct_ancestry.csv.gz");
magicancestry = readmagicancestry("outstem_magicreconstruct_ancestry.csv.gz");
truefgl = formmagicgeno(string(dataid,"_magicsimulate_truefgl.csv.gz"),pedinfo; isfounderinbred);
acc = magicaccuracy!(truefgl, magicancestry; isfounderinbred)
println(acc)

plotcondprob(magicancestry;
    truefgl, 
    offspring=1,
    probtype="haploprob"
)




# clear up
# cd(@__DIR__)
# dataid = "sim"
# rm.(filter(x->occursin(dataid, x), readdir()))

