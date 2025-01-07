using MagicScan
cd(@__DIR__)
pwd()


using CSV, DataFrames
using StatsModels
using StatsBase

ancestryfile = outstem*"_magicreconstruct_ancestry.csv.gz"
phenofile = outstem*"_magicsimulate_pheno.csv"
phenodf = CSV.read(phenofile,DataFrame)
f = @formula(phenotype ~ 1)
f2 = apply_schema(f, schema(f, phenodf))
yresp,xpred = modelcols(f2,phenodf)

magicancestry = readmagicancestry(ancestryfile);
noff = size(magicancestry.magicped.offspringinfo,1)
haploprob = MagicScan.get_haploprob(magicancestry, 1:noff)
loglikels,lodls,log10pls = MagicScan.scan_lm(yresp,xpred,haploprob);

nperm = 100
res = zeros(nperm,3)
ny = length(yresp)
i = 1
y = sample(yresp,ny,replace=false)
loglikels,lodls,log10pls = MagicScan.scan_lm(y,xpred,haploprob);
lodls
