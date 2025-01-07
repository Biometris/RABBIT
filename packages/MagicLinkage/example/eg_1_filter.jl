
using Revise
using MagicFilter
cd(@__DIR__)
pwd()

isfounderinbred = true 

# pre-step filter markers
dataid = "sim"
genofile=string(dataid,"_magicsimulate_geno.vcf.gz")
pedfile = string(dataid,"_magicsimulate_ped.csv")
outstem = dataid
@time magicfilter(genofile,pedfile;
    isfounderinbred,
    snp_mono2miss = true,    
    snp_missfilter = (f,o) -> o < 0.8,
    snp_minmaf = 0.1,    
    outstem 
)
