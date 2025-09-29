
# using Distributed
# addprocs(5) # add worker processes on local machine
# @info string("nprocs=", nprocs())
# @everywhere  using MagicFilter

using Revise
using MagicBase, MagicFilter
cd(@__DIR__)
pwd()
using CSV, DataFrames

genofile = "sim_magicsimulate_geno.vcf.gz"
pedfile = "sim_magicsimulate_ped.csv"
magicgeno = formmagicgeno(genofile,pedfile);

@time MagicFilter.filter_marker!(magicgeno;
    mono2miss = true,    
);

# cd(@__DIR__)
# outid = "sim"
# rm.(filter(x->occursin(outid,x), readdir()))

fls = [0.005, 0.0359, 0.0669, 0.0978, 0.1288, 0.3762, 0.4072, 0.4381, 0.4691, 0.5, 0.5309, 0.5619, 0.5928, 0.6237, 0.8712, 0.9022, 0.9331, 0.9641, 0.995]

fls =  [0.005, 0.5, 0.995]

n1 = 0
n2 = 6
MagicFilter.cal_post_prob(n1,n2,fls)