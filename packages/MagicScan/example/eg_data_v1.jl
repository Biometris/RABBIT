using MagicBase, MagicSimulate, MagicReconstruct
cd(@__DIR__)


# simulate founder haplotypes
dataid = "sim"
fhaplofile = dataid*"_fhaplo.vcf.gz"
ncluster = 5
@time simfhaplo(
    nsnp = 100*ncluster,
    nparent = 4,
    chrlen = 100*ones(ncluster),
    outfile = fhaplofile
)

# simulate genotypic/phenotypic data
magicsimulate(fhaplofile,"4ril-self6";
    popsize = 300,
    ispheno = true,
    pheno_nqtl = 1,
    pheno_h2 = 0.2,
    outstem = dataid
)

truepheno = MagicBase.readmultitable(dataid*"_magicsimulate_truepheno.csv")
println("true QTL:\n",truepheno["map_qtl"][!,1:5])

# haplotypoe reconstruct
genofile = dataid*"_magicsimulate_geno.vcf.gz"
magicreconstruct(genofile,"4ril-self6"; 
    model = "depmodel",
    outstem = dataid
)

#  clear up
rm(fhaplofile)
rm.(filter(x->occursin(".log",x), readdir()))
rm.(filter(x->occursin("geno",x), readdir()))
rm.(filter(x->occursin("_ped",x), readdir()))
rm.(filter(x->occursin("fgl",x), readdir()))

