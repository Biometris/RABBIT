using Revise
using MagicBase
cd(@__DIR__)
pwd()

dataid = "outstem"
vcffile = dataid*"_geno.vcf"
pedfile = dataid*"_ped.csv"
@time magicgeno = formmagicgeno(vcffile,pedfile);

@time parsevcf(vcffile, pedfile; 
    formatpriority=["GT","AD"],
    outstem = dataid*"_afterparse"
)


rm.(filter!(x->occursin("outstem",x),readdir()))


genofile = "f2_F0173_F0606_magicimpute_geno.vcf"