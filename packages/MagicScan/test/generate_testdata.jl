using MagicBase, MagicSimulate, MagicReconstruct
cd(@__DIR__)


# simulate founder haplotypes
dataid = "sim"
fhaplofile = dataid*"_fhaplo.csv.gz"
ncluster = 2
@time simfhaplo(
    nsnp = 150*ncluster,
    nparent = 4,
    chrlen = 100*ones(ncluster),
    outfile = fhaplofile
)

# simulate genotypic/phenotypic data
magicsimulate(fhaplofile,"4ril-self6";
    popsize = 200,
    ispheno = true,
    pheno_nqtl = 1,
    pheno_h2 = 0.5,
    outstem = dataid
)

# haplotypoe reconstruct
genofile = dataid*"_magicsimulate_geno.vcf.gz"
magicreconstruct(genofile,"4ril-self6";
    model = "depmodel",
    outstem = dataid
)

#  clean up
rm(fhaplofile)
rm.(filter(x->occursin(".log",x), readdir()))
rm.(filter(x->occursin("geno",x), readdir()))
rm.(filter(x->occursin("_ped",x), readdir()))
rm.(filter(x->occursin("fgl",x), readdir()))
rm.(filter(x->occursin("_plots.tar",x), readdir()))
rm.(filter(x->occursin("_recom.",x), readdir()))
