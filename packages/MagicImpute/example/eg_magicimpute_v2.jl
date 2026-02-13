

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
    magicgeno.foundergeno[chr][:,1:2] .= missingcode
end
@time magicimpute!(magicimpute; 
    isfounderinbred,           
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


genofile2 = outstem*"_geno.vcf.gz"
@time magicancestry=magicreconstruct(genofile,pedfile;
    # model = "depmodel",    
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




# clear up
# cd(@__DIR__)
# dataid = "sim"
# rm.(filter(x->occursin(dataid, x), readdir()))

