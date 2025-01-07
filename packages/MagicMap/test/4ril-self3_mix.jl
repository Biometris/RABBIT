using MagicBase
using MagicMap
using Test
cd(@__DIR__)

# construct map
designcode = "4ril-self3"
dataid = designcode*"_mix"
genofile = dataid*"_magicsimulate_geno.vcf.gz"
outstem=dataid
ncluster = 2
outfiles = magicmap(genofile, designcode; 
    ncluster,
    outstem,
    verbose=false
)

truemapfile = string(dataid,"_magicsimulate_truegeno.csv.gz")
estmapfile = outstem*"_magicmap_construct_map.csv.gz"
acc = mapaccuracy(truemapfile,estmapfile; isgroupacc=true)

# println(acc)
@test sum(abs.(acc.orderacc))/length(acc.orderacc) > 0.6

#  clean up
rm.(filter(x->occursin(outstem*"_magicmap",x), readdir()))
