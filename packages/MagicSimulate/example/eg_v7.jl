using Revise
using MagicSimulate
using MagicBase
using Distributions
cd(@__DIR__)
pwd()



@time simfhaplo(nsnp=1000,nparent=20, 
    isfounderinbred = false,
    outfile="sim_fhaplo.vcf"
)

# matescheme = MateScheme(2,["Pairing","DH"],[1,1])
# matescheme = MateScheme(4,["Pairing","DH"],[2,1])
# matescheme = MateScheme(4,["Pairing","Selfing"],[2,6])
# matescheme = MateScheme(8,["Pairing","Selfing"],[2,4])
# matescheme = MateScheme(8,["Pairing","Sibling"],[2,10])
# matescheme = MateScheme(19,["FullDiallel", "RM1_E", "Selfing"],[1,4,6])

pedinfo = "4star-self1"

@time magicsimulate("sim_fhaplo.vcf",pedinfo;
    popsize=200,    
    isfounderinbred = false,
    outstem="sim"
)

rm.(filter(x->occursin("sim",x), readdir()))


