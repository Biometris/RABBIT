
using Revise
using MagicBase, MagicSimulate
using Distributions
cd(@__DIR__)
pwd()

fhaplofile = "sim_fhaplo.csv.gz"
nchr = 10
simfhaplo(nsnp=200*nchr,nparent=4,
    isfounderinbred=true,
    chrlen = 100*ones(nchr),
    outfile= fhaplofile, 
)


nparent = 4
popsize =  200
epsf = epso = 0.01
missf = misso = 0.2
pedinfo = string(nparent, "ril-self6")
@time res = [begin 
    simres = magicsimulate(fhaplofile,pedinfo;
        foundererror = Beta(2,1/epsf-2),
        offspringerror = Beta(2,1/epso-2),
        foundermiss = Beta(2,1/missf-2),
        offspringmiss = Beta(2,1/misso-2),      
        popsize, 
        select_nqtl = 100, 
        select_dom,
        select_prop = 0.1,    
        outstem=nothing,
    )
    fhetero = mean([mean(allunique.(i)) for i in simres.truegeno.offspringgeno])
    nonibd = mean([mean(allunique.(i)) for i in simres.magicfgl.offspringgeno])
    [select_dom,fhetero,nonibd]
end for select_dom in [0,1,3,5,7,9,11]]
# for select_prop in [1.0,0.75,0.5,0.2,0.1,0.05]]

using Plots
x = first.(res)
y = [i[2] for i in res]
z = [i[3] for i in res]
plot(scatter(x,y),scatter(x,z),layout=(2,1))