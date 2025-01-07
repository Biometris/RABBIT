using MagicBase, MagicSimulate, MagicFilter
using CSV, DataFrames
cd(@__DIR__)


# simulate founder haplotypes
fhaplofile = "fhaplo.vcf.gz"
ncluster = 2
nsnpchr = 200
simfhaplo(nsnp=ncluster*nsnpchr, nparent=6,
    chrlen=100*ones(ncluster),
    outfile=fhaplofile
)

magicped = MagicBase.generate_magicped(;
    designcodes = ["bc3-self2","P2/P3=>DH","P4/3/P2/P3//P5/P6=>3","2ril-self6"],
    founders = ["P1||P2","NA", "NA", "P5||P6"],
    subpopsizes =50*ones(Int, 4));
pedfile = "multipop_ped.csv"
savemagicped(pedfile, magicped)
# plotmagicped(magicped)

outstem = "multipop"
magicsimulate(fhaplofile,pedfile; outstem)

magicfilter(outstem*"_magicsimulate_geno.vcf.gz", 
    outstem*"_magicsimulate_ped.csv";  
    outstem
)
    
#  clean up
rm(fhaplofile)
rm.(filter(x->occursin(".log",x), readdir()))
rm.(filter(x->occursin("_ped.png",x), readdir()))
rm.(filter(x->occursin("_truecontfgl",x), readdir()))
rm.(filter(x->occursin("_truefgl",x), readdir()))
rm.(filter(x->occursin("_magicfilter_snp",x), readdir()))
rm.(filter(x->occursin("_magicfilter_ind",x), readdir()))
rm.(filter(x->occursin("_magicfilter_founder",x), readdir()))
rm.(filter(x->occursin("_magicfilter_subpop",x), readdir()))