
using MagicBase, MagicScan
using CSV,DataFrames
cd(@__DIR__)
pwd()

dataid = "sim"
ancestryfile = dataid*"_magicreconstruct_ancestry.csv.gz"
phenofile = dataid*"_magicsimulate_pheno.csv"
outstem=dataid
magicscan(ancestryfile,phenofile; 
    thresholds = nothing,    
    equation = @formula(phenotype ~ 1), 
    outstem,    
)
filter(x->occursin(outstem,x), readdir())

truepheno = MagicBase.readmultitable(dataid*"_magicsimulate_truepheno.csv");
println("trueqtl: ", truepheno["map_qtl"][!,1:5])

CSV.read(outstem*"_magicscan_peak.csv",DataFrame)

# rm.(filter!(x->occursin(dataid,x),readdir()))

# rm.(filter!(x->occursin(outstem,x),readdir()))
