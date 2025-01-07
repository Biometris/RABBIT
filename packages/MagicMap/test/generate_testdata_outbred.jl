
using MagicBase, MagicSimulate, MagicFilter
cd(@__DIR__)

# simulate data
designcode = "4star-self1"
fhaplofile = "fhaplo.vcf"
ncluster = 2
nsnpchr = 200
simfhaplo(nsnp=ncluster*nsnpchr, nparent=8,
    isfounderinbred = false,
    chrlen=100*ones(ncluster),
    outfile=fhaplofile
)
dataid = designcode*"_outbred"
magicsimulate(fhaplofile,designcode;
    popsize = 150,
    isfounderinbred = false,
    outstem = dataid,
)

genofile=string(dataid,"_magicsimulate_geno.vcf.gz")
pedfile = string(dataid,"_magicsimulate_ped.csv")
@time magicfilter(genofile,pedfile;
    isfounderinbred = false,
    snp_missfilter = (f,o) -> o < 0.8,
    snp_minmaf = 0.1,    
    outstem = dataid, 
)

rm(fhaplofile)
rm(dataid*"_magicsimulate_geno.vcf.gz")
rm.(filter(x->occursin(".log",x), readdir()))
rm.(filter(x->occursin("_ped.png",x), readdir()))
rm.(filter(x->occursin("_truecontfgl",x), readdir()))
rm.(filter(x->occursin("_truefgl",x), readdir()))
files = setdiff(readdir(),[dataid*"_magicfilter"*i for i in ["_geno.vcf.gz","_ped.csv"]])
filter!(x->occursin("_magicfilter",x),files)
rm.(files)
