using Revise
using MagicBase, MagicSimulate, MagicReconstruct
using  Distributions
cd(@__DIR__)
pwd()


for isinbred in [true, false]
    fhaplofile = isinbred ? "fhaplo_inbred.csv" : "fhaplo_outbred.csv"
    ncluster = 5
    simfhaplo(nsnp=100*ncluster,
        nparent= 6,
        isfounderinbred = isinbred,
        chrlen=100*ones(ncluster),
        outfile=fhaplofile
    )
end


designcode = "4ril-self3"
magicsimulate("fhaplo_inbred.csv",designcode;
    popsize = 5,
    outstem="simarray"    
)

magicsimulate("fhaplo_inbred.csv",designcode;
    popsize = 5,
    outstem="simmix",
    seqfrac = 0.5,
    seqdepth = Gamma(2,10),
)

magicped = MagicBase.generate_magicped(;
    designcodes = ["bc3-self2","P2/P3=>DH","P4/3/P2/P3//P5/P6=>3","2ril-self6"],
    founders = ["P1||P2","NA", "NA", "P5||P6"],
    subpopsizes =2*ones(Int, 4));
pedfile = "multipop_ped.csv"
savemagicped(pedfile, magicped)
magicsimulate("fhaplo_inbred.csv", pedfile;    
    outstem = "multipop",
)


magicreconstruct("simarray_magicsimulate_geno.vcf.gz",designcode;
    chrsubset=[2,4],
    outstem="simarray",
)

magicreconstruct("simmix_magicsimulate_geno.vcf.gz",designcode;
    chrsubset=[1,5],
    outstem="simmix",
)

magicreconstruct("multipop_magicsimulate_geno.vcf.gz",pedfile;    
    chrsubset=[3],
    outstem="multipop",
)

rm.(filter(x->occursin(".log",x), readdir()))
rm.(filter(x->occursin(".png",x), readdir()))
rm.(filter(x->occursin("_plots.tar",x), readdir()))
rm.(filter(x->occursin("_recom.csv",x), readdir()))
