using Revise
using MagicSimulate
using MagicBase
using Distributions
cd(@__DIR__)
pwd()

@time simfhaplo(nsnp=1000,nparent=8, outfile="sim_fhaplo.csv")

@time magicsimulate("sim_fhaplo.csv","4ril-self3";
    popsize=200,    
    seqfrac = 0.5,
    outstem="sim"
)



@time obsgeno, truegeno, magicfgl, contfgl, phenores = magicsimulate(
    "sim_fhaplo.csv","8ril-self4";
    popsize=200,        
    # error_randallele = 0,
    outstem=nothing
);

# rm.(filter(x->occursin("simarray",x), readdir()))

chr = 1
tt = join.(truegeno.offspringgeno[chr],"")
tt[tt .== "21"] .= "12"

obs = obsgeno.offspringgeno[chr]
b =  obs .!== tt
b2 = obs .!== "NN"
b .*= b2
a = [tt[b] obs[b]]
unique(a[:,1])

ls = collect(eachrow(a))
sum([i == ["11","22"] for i in ls])
sum([i == ["11","12"] for i in ls])
sum([i == ["22","11"] for i in ls])
sum([i == ["22","12"] for i in ls])
sum([i == ["12","11"] for i in ls])
sum([i == ["12","22"] for i in ls])
