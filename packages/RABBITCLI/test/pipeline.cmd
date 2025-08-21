:: set WORKDIR and RABBITDIR 
set "WORKDIR=%~dp0"
cd "%WORKDIR%"
set "RABBITDIR=%WORKDIR%..\src\"
set OUTSTEM=sim

julia "%RABBITDIR%rabbit_simfhaplo.jl" --nsnp 500 --nparent 3 ^
  --chrlen "[100,100]" ^
  -o %OUTSTEM%_fhaplo.vcf.gz

julia "%RABBITDIR%rabbit_generate_magicped.jl" --designcodes "[P1/P2=>DH, 2ril-self1]" ^
  --founders "[NA,P1||P3]" ^
  --subpopsizes "[100,100]" ^
  -o %OUTSTEM%

julia "%RABBITDIR%rabbit_magicsimulate.jl" -g %OUTSTEM%_fhaplo.vcf.gz -p %OUTSTEM%_ped.csv ^
  --popsize 200 --seqfrac 1.0  ^
  --foundermiss Beta(1,9) --offspringmiss Beta(1,9) ^
  --foundererror Beta(1,19) --offspringerror Beta(1,19) ^
  --seqdepth Gamma(2,5) ^
  --ispheno true --pheno_nqtl 1 --pheno_h2 0.5 ^
  -o %OUTSTEM%

julia "%RABBITDIR%rabbit_magicfilter.jl" -g %OUTSTEM%_magicsimulate_geno.vcf.gz -p %OUTSTEM%_magicsimulate_ped.csv ^
  --maxomiss 0.7 --minmaf 0.05 --mono2miss true ^
  -o %OUTSTEM%

julia "%RABBITDIR%rabbit_magiccall.jl" -g %OUTSTEM%_magicfilter_geno.vcf.gz -p %OUTSTEM%_magicfilter_ped.csv ^
  --nworker 2 -o %OUTSTEM%

julia "%RABBITDIR%rabbit_magicmap.jl" -g %OUTSTEM%_magiccall_geno.vcf.gz -p %OUTSTEM%_magicfilter_ped.csv ^
  --ncluster 2 ^
  --nworker 2 -o %OUTSTEM%

julia "%RABBITDIR%rabbit_magicmask_impute.jl" -g %OUTSTEM%_magicfilter_geno.vcf.gz -p %OUTSTEM%_magicfilter_ped.csv ^
  --mapfile %OUTSTEM%_magicmap_construct_map.csv.gz --isordermarker false ^
  --nworker 2 -o %OUTSTEM%

julia "%RABBITDIR%rabbit_magicreconstruct.jl" -g %OUTSTEM%_magicimpute_geno.vcf.gz -p %OUTSTEM%_magicfilter_ped.csv ^
  --nworker 2 -o %OUTSTEM%

julia "%RABBITDIR%rabbit_magicscan.jl" -g %OUTSTEM%_magicreconstruct_ancestry.csv.gz -p %OUTSTEM%_magicsimulate_pheno.csv ^
  -o %OUTSTEM%


