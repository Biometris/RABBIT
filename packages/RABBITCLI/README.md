# RABBITCLI

RABBITCLI consists of the wrappers of RABBIT functions, such that these functions can be run as command lines in a command shell.

## Get started

1. Download and install Julia available at https://julialang.org/
2. Download RABBIT repository
3. In a command shell, change into work directory including input files, and for example, run command line:

```
path/to/bin/julia.exe path/to/RABBITCLI/src/rabbit_funcname.jl -g geno.csv -p ped.csv
```

## Help

```
path/to/bin/julia.exe path/to/RABBITCLI/src/rabbit_funcname.jl --help
```

## Main function files

* `rabbit_simfhaplo.jl` simulates founder haplotypes
* `rabbit_generate_magicped.jl` generates pedfile. 
* `rabbit_magicsimulate.jl` simulates offspring genotypic data
* `rabbit_parsebreedped.jl` parses BASF-format pedigree file
* `rabbit_vcffilter.jl` filters vcf genofile
* `rabbit_magicparse.jl` parses genofile and pedfile
* `rabbit_magicfilter.jl` filters markers and individuals
* `rabbit_magiccall.jl` single marker genotype calling
* `rabbit_magicmap.jl` map construction
* `rabbit_magicmask.jl` masks genotypes for imputation evaluation
* `rabbit_magicimpute.jl` genotype imputation
* `rabbit_imputeaccuracy.jl` calculates imputation accuracy
* `rabbit_magicmask_impute.jl` combines magicmask, magicimpute, and imputeaccuracy. 
* `rabbit_magicreconstruct.jl` haplotype reconstruction
* `rabbit_magicscan.jl` QTL scan using a linear model
