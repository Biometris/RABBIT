
using MagicBase, MagicSimulate, MagicImpute
cd(@__DIR__)

# simulate data
designcode = "4ril-self3"
fhaplofile = "fhaplo.vcf"
ncluster = 1
nsnpchr = 200
simfhaplo(nsnp=ncluster*nsnpchr, nparent=8,
    isfounderinbred = false,
    chrlen=100*ones(ncluster),
    outfile=fhaplofile
)
magicsimulate(fhaplofile,designcode;
    popsize = 50,
    isfounderinbred = false,
    outstem = designcode*"_outbred",
)

dataid = designcode*"_outbred"
genofile = dataid*"_magicsimulate_geno.vcf.gz"
magicimpute(genofile,designcode;
    isfounderinbred = false, 
    outstem = dataid,
    # isdelmarker = true, 
    # iscorrectfounder = true, 
    # isinfererror = false,
    # verbose=false
)


#  clean up
rm.(filter(x->occursin(dataid*"_magicimpute_founder",x), readdir()))
rm.(filter(x->occursin("genofreq",x) || occursin("_map.csv",x), readdir()))
rm.(filter(x->occursin("inferred_error",x), readdir()))
rm(fhaplofile)
rm.(filter(x->occursin(".log",x), readdir()))
rm.(filter(x->occursin("_ped",x) && occursin(designcode,x), readdir()))
rm.(filter(x->occursin("_truecontfgl",x) || occursin("_truegeno",x), readdir()))
