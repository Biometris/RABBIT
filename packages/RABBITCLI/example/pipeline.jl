using Distributed
nprocs() < 5 && addprocs(5-nprocs()) # 4 chromosomes
@info string("nworkers=", nworkers())
@everywhere using MagicFilter, MagicCall, MagicMap, MagicImpute, MagicReconstruct, MagicScan

using Revise
using MagicFilter, MagicCall, MagicMap, MagicImpute, MagicReconstruct
using MagicBase, MagicSimulate, MagicScan
cd(@__DIR__)
outstem = "example"
# TEST scenarios: istart, isfounderinbred, seqfrac =0/1

isstar = false
isfounderinbred = !isstar
# isfounderinbred = true 
# S1_1 Simulate founder genotypic data
fhaplofile = outstem*"_fhaplo.vcf"    
simfhaplo(;
    isfounderinbred,
    nsnp = 400, 
    nparent = isstar ? 3 : 4,
    chrlen = 100*ones(4),
    outfile = fhaplofile
)

# S1_2 Simulate pedigree information
if isstar
    magicped = formmagicped("3star-self1",200)
else
    magicped = generate_magicped(;
    designcodes=["P1/P2=>DH", "2ril-self1", "P4/3/P4//P2/P3=>1"],
    founders = ["NA","P1||P3","NA"],
    subpopsizes=100*ones(3)
    )
end
savemagicped(outstem*"_ped.csv", magicped)
plotmagicped(magicped)

# S1_3 Simulate offspring genotypic data
using Distributions
pedfile = outstem*"_ped.csv"
magicsimulate(fhaplofile,pedfile;    
    isfounderinbred,
    seqfrac = 1.0,
    seqdepth = Gamma(2,5),    
    foundermiss = Beta(1,9),
    offspringmiss = Beta(1,9),
    foundererror = Beta(1,19),
    offspringerror = Beta(1,19),    
    allelebalancemean = Beta(5,5),
    allelebalancedisperse = Exponential(0.05),    
    ispheno = true,
    pheno_nqtl=1,
    pheno_h2= 0.5,
    nplot_subpop = 10, 
    outstem,
)   


# S2 Data filtering
genofile = outstem*"_magicsimulate_geno.vcf.gz"
pedfile = outstem*"_magicsimulate_ped.csv"
magicfilter(genofile,pedfile;        
    isfounderinbred,
    minmaf = 0.05,
    missfilter = (f,o) -> o <= 0.7,     
    offspring_maxmiss = 0.95,
    # isfilterdupe = true,  
    outstem
);


# S3 genotype calling
genofile = outstem*"_magicfilter_geno.vcf.gz"
pedfile = outstem*"_magicfilter_ped.csv"
magiccall(genofile,pedfile;   
    isfounderinbred,               
    isparallel = false,     
    outstem 
)
truefile = outstem*"_magicsimulate_truegeno.csv.gz"
calledgenofile = outstem*"_magiccall_geno.vcf.gz"
acc = magicaccuracy(truefile,calledgenofile,pedfile;isfounderinbred)
println(acc)
plotmarkererror(calledgenofile)


# S4 map construction 
# genofile = outstem*"_magicfilter_geno.vcf.gz"
genofile = outstem*"_magiccall_geno.vcf.gz"
pedfile = outstem*"_magicfilter_ped.csv"
magicmap(genofile,pedfile;        
    isfounderinbred,               
    ncluster = 4,              
    outstem
)

# S5 genotype imputation
genofile = outstem*"_magicfilter_geno.vcf.gz"
pedfile = outstem*"_magicfilter_ped.csv"
magicgeno = magicimpute(genofile,pedfile;
    isfounderinbred,               
    mapfile = outstem*"_magicmap_construct_map.csv.gz",                     
    outstem         
); 

# S6 haplotype reconstruct
genofile = outstem*"_magicimpute_geno.vcf.gz"
pedfile = outstem*"_magicfilter_ped.csv"
magicancestry = magicreconstruct(genofile,pedfile;         
    isfounderinbred,               
    nplot_subpop = 1, 
    outstem     
);

truefgl = formmagicgeno(outstem*"_magicsimulate_truefgl.csv.gz",pedfile;isfounderinbred);
acc = magicaccuracy!(truefgl, magicancestry; isfounderinbred) 
println(acc)
fig = plotcondprob(magicancestry; truefgl,
    probtype="genoprob", 
    offspring=1,
)
display(fig)

# S7 QTL scan
using StatsModels # required for @formula
ancestryfile = outstem*"_magicreconstruct_ancestry.csv.gz"
phenofile = outstem*"_magicsimulate_pheno.csv"
peak = magicscan(ancestryfile,phenofile;
    equation = @formula(phenotype ~ 1 + population), 
    outstem     
)
println("Profile peak: \n", peak) 
truepheno = MagicBase.readmultitable(outstem*"_magicsimulate_truepheno.csv");
println("trueqtl: ", truepheno["map_qtl"])

# clean up
# cd(@__DIR__)
# rm.(filter(x->occursin("example_", x),readdir()))
# cd(@__DIR__)
# rm.(filter(x->occursin("sim_", x),readdir()))


