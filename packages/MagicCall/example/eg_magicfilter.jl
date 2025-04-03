
using Revise
using MagicFilter
cd(@__DIR__)
pwd()


dataid = "sim"
genofile=string(dataid,"_magicsimulate_geno.vcf.gz")
pedfile = string(dataid,"_magicsimulate_ped.csv")
outstem = dataid
@time magicfilter(genofile,pedfile;
    isfounderinbred = true,
    mono2miss = true,    
    missfilter = (f,o) -> o < 0.9,
    minmaf = 0.05,    
    outstem 
)