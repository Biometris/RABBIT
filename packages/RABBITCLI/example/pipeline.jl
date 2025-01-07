using Distributed
nprocs() < 5 && addprocs(5-nprocs()) # 4 chromosomes
@info string("nworkers=", nworkers())
@everywhere using MagicFilter, MagicCall, MagicMap, MagicImpute, MagicReconstruct

using MagicFilter, MagicCall, MagicMap, MagicImpute, MagicReconstruct
using MagicBase, MagicSimulate, MagicScan
cd(@__DIR__)
outstem = "example"
isstar = true

# S1_1 Simulate founder genotypic data
fhaplofile = outstem*"_fhaplo.vcf.gz"    
simfhaplo(;
    isfounderinbred = true,
    nsnp = 400, 
    nparent = isstar ? 5 : 4,
    chrlen = 100*ones(4),
    outfile = fhaplofile
)

# S1_2 Simulate pedigree information
if isstar
    magicped = formmagicped("5star-self1",400)
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
    snp_minmaf = 0.05,
    snp_missfilter = (f,o) -> o <= 0.7,     
    offspring_maxmiss = 0.95,
    isfilterdupe = true,  
    outstem
);


# S3 genotype calling
genofile = outstem*"_magicfilter_geno.vcf.gz"
pedfile = outstem*"_magicfilter_ped.csv"
magiccall(genofile,pedfile;               
    outstem 
)

truefile = outstem*"_magicsimulate_truegeno.csv.gz"
calledgenofile = outstem*"_magiccall_geno.vcf.gz"
acc = magicaccuracy(truefile,calledgenofile,pedfile)
println(acc)
plotmarkererror(calledgenofile)

# S4 map construction 
genofile = outstem*"_magiccall_geno.vcf.gz"
pedfile = outstem*"_magicfilter_ped.csv"
magicmap(genofile,pedfile;        
    minncluster = 2, 
    maxncluster = 10,                      
    outstem
)


# S5 genotype imputation
genofile = outstem*"_magicfilter_geno.vcf.gz"
pedfile = outstem*"_magicfilter_ped.csv"
magicgeno = magicimpute(genofile,pedfile;
    mapfile = outstem*"_magicmap_construct_map.csv.gz",                     
    outstem         
); 

# S6 haplotype reconstruct
genofile = outstem*"_magicimpute_geno.vcf.gz"
pedfile = outstem*"_magicfilter_ped.csv"
magicancestry = magicreconstruct(genofile,pedfile;         
    nplot_subpop = 1, 
    # formatpriority = ["GT"],    
    outstem     
);

truefgl = formmagicgeno(outstem*"_magicsimulate_truefgl.csv.gz",pedfile);
acc = magicaccuracy!(truefgl, magicancestry) 
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