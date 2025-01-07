using MagicBase, MagicScan
using Test
cd(@__DIR__)
pwd()

dataid = "sim"
ancestryfile = dataid*"_magicreconstruct_ancestry.csv.gz"
@test isa(readmagicancestry(ancestryfile),MagicAncestry)
phenofile = dataid*"_magicsimulate_pheno.csv"
@time peak = magicscan(ancestryfile,phenofile;
    thresholds = nothing, nperm = 100, 
    verbose = false,
)

map_qtl = MagicBase.readmultitable(dataid*"_magicsimulate_truepheno.csv")["map_qtl"]
qtl_true = (string(map_qtl[1,:linkagegroup]), map_qtl[1,:poscm])
qtl_est_int = parse.(Float64,split(peak[6,2],"~"))
qtl_est_chr = peak[2,2]

@test qtl_est_chr == qtl_true[1]
@test qtl_est_int[1]-1.0 <= qtl_true[2] <= qtl_est_int[2]+1.0

# clearnup
rm.(filter!(x->occursin("outstem",x),readdir()))
