
using Revise
using MagicBase, MagicSimulate
using Distributions
cd(@__DIR__)
pwd()



@time obsgeno, truegeno, magicfgl, contfgl,phenores= magicsimulate(fhaplofile,"4ril-self2";
    outstem=nothing,
    select_nqtl,select_dom,select_prop,
    ispheno =true
);


using MagicSimulate
using Distributions
cd(@__DIR__)
pwd()

nparent = 4 
popsize =  200
pedinfo = string(nparent, "ril-self4")
println(outstem)
res = magicsimulate(fhaplofile,pedinfo;
    foundererror = Beta(2,1/epsf-2),
    offspringerror = Beta(2,1/epso-2),
    foundermiss = Beta(2,1/missf-2),
    offspringmiss = Beta(2,1/misso-2),      
    popsize, 
    select_prop,
    select_nqtl = 10, 
    outstem, workdir
)