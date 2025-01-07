
using MagicSimulate
using MagicBase
using Distributions
using Test
cd(@__DIR__)


# fhaplofile = "fhaplo_inbred.vcf.gz"
# outfile = simfhaplo(nsnp=200,nparent=4, chrlen = 100*ones(2), outfile=fhaplofile)
# @test relpath(outfile,fhaplofile) == "."
# outfile = MagicBase.tocsvgeno(fhaplofile,4)
# @test relpath(outfile,replace(fhaplofile,".vcf"=>".csv")) == "."

@time @testset "MagicSimulate" begin
    @testset "pedfile" begin
        @test isa(magicsimulate("fhaplo_inbred.csv.gz","ped.csv";
            popsize = 5,
            outstem="simarray",
            verbose=false,
        ), Vector{String})
        @test isa(magicsimulate("fhaplo_inbred.vcf.gz","ped.csv";
            popsize = 5,
            seqfrac = 0.5,
            outstem="simgbs",
            verbose=false,
        ), Vector{String})
    end
    @testset "simfhgaplo" begin
        @test isfile(simfhaplo(nsnp=100,nparent=4,
            isfounderinbred=true,
            outfile="simarray_fhaplo.csv",            
        ))
        @test isfile(simfhaplo(nsnp=100,nparent=4,
            isfounderinbred=false,
            outfile="simarray_fphased.csv",            
        ))
    end
    @testset "4ril-self3" begin
        @test isa(magicsimulate("fhaplo_inbred.vcf.gz","4ril-self3";
            popsize = 5,
            foundererror = Beta(1,99),
            offspringerror = Beta(1,99),
            foundermiss = Beta(1,9),
            offspringmiss = Beta(1,9),
            isobligate = true,
            interference = 5,
            outstem=nothing,
            verbose=false,
        ), NamedTuple)
        @test isa(magicsimulate("fhaplo_inbred.csv.gz","4ril-self3";
            popsize = 5,
            outstem="simarray",
            verbose=false,
        ), Vector{String})
        @test isa(magicsimulate("fhaplo_inbred.vcf.gz","4ril-self3";
            popsize = 5,
            seqfrac = 0.5,
            seqerror = Beta(1,199),
            seqdepth = Gamma(2,10),
            isobligate = true,
            interference = 5,
            outstem=nothing,
            verbose=false,
        ), NamedTuple)
        @test isa(magicsimulate("fhaplo_inbred.csv.gz","4ril-self3";
            popsize = 5,
            seqfrac = 0.5,
            outstem="simgbs",
            verbose=false,
        ), Vector{String})
    end
    @testset "output" begin
        # obsgeno, truegeno, magicfgl, contfgl,phenores = res
        res = magicsimulate("fhaplo_inbred.csv.gz","4ril-self3";
            popsize = 5,
            ispheno = true,
            outstem=nothing,
            verbose=false,
        )
        for i in res
            isnothing(i) && continue
            isa(i,NamedTuple) && continue
            isa(i,MagicPed) && continue
            @test isa(savegenodata("test.csv",i),String)
            @test isa(savemagicgeno("test.csv",i),String)
            @test isa(readmagicgeno("test.csv"),MagicGeno)
        end
        res = magicsimulate("fhaplo_inbred.vcf.gz","4ril-self3";
            seqfrac = 0.5,
            ispheno = true,
            popsize = 5,
            outstem=nothing,
            verbose=false,
        )
        for i in res
            isnothing(i) && continue
            isa(i,NamedTuple) && continue
            isa(i,MagicPed) && continue
            @test isa(savegenodata("test.csv",i),String)
            @test isa(savemagicgeno("test.csv",i),String)
            @test isa(readmagicgeno("test.csv"),MagicGeno)
        end
    end
end
rm("test.csv")
# rm("fhaplo_inbred.vcf.gz")
rm.(filter(x->occursin("simarray",x), readdir()))
rm.(filter(x->occursin("simgbs",x), readdir()))
