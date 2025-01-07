
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
    snp_mono2miss = true,    
    snp_missfilter = (f,o) -> o < 0.9,
    snp_minmaf = 0.05,    
    outstem 
)