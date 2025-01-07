# using Revise
using MagicBase
using Test
cd(@__DIR__)

verbose = false

@time @testset "MagicBase" begin
    verbose && println("1")
    @testset "MagicPed" begin
        codels = ["8ril-self4","8ril-sib10","6star-self6","bc2-self3"]
        popsize = 10
        magicpedfile = "magicped.csv"
        for designcode = codels
            @test isa(formmagicped(designcode,popsize),MagicPed)
            magicped = formmagicped(designcode,popsize)
            @test isa(plotmagicped(magicped),Any)
            @test isfile(savemagicped(magicpedfile,magicped))
            @test isa(readmagicped(magicpedfile), MagicPed)
        end
        rm(magicpedfile)
    end
    verbose && println("2")
    @testset "MagicGeno" begin
        # form/save/read/ magicgeno
        for dataid = ["simarray_magicsimulate","simmix_magicsimulate","multipop_magicsimulate"]
            genofile = dataid*"_geno.vcf.gz"
            @test isa(formmagicgeno(genofile,4),MagicGeno)
            @test isa(formmagicgeno(genofile,dataid*"_ped.csv"),MagicGeno)
            magicgeno = formmagicgeno(genofile,dataid*"_ped.csv")
            @test isfile(savemagicgeno("magicgeno.csv",magicgeno))
            @test isa(readmagicgeno("magicgeno.csv"),MagicGeno)
            @test isa(MagicBase.submagicgeno!(magicgeno,
                chrsubset=[2,4],snpsubset=1:3:1000),MagicGeno)
            @test isfile(savegenodata("geno.csv",magicgeno))
            @test isa(formmagicgeno("geno.csv",dataid*"_ped.csv"), MagicGeno)
            @test isa(MagicBase.submagicgeno(magicgeno,
                chrsubset=[2,4],snpsubset=1:3:1000),MagicGeno)
            @test isa(MagicBase.rawgenoprob!(magicgeno), MagicGeno)
            @test isa(MagicBase.rawgenocall!(magicgeno),MagicGeno)
        end
        rm("magicgeno.csv")
        rm("geno.csv")
    end
    verbose && println("3")
    @testset "formfhaplo" begin
        for isinbred = [false,true]
            fhaplofile = isinbred ? "fhaplo_inbred.csv" : "fhaplo_outbred.csv"
            @test isa(MagicBase.formfhaplo(fhaplofile,
                isfounderinbred=isinbred), MagicGeno)
            fhaplo = MagicBase.formfhaplo(fhaplofile,
                isfounderinbred=isinbred)
            if isinbred
                fhaplo.foundergeno[1][1,2] = "1"
                fhaplo.foundergeno[1][2,3] = "N"
                fhaplo.foundergeno[1][3,4] = "N"
            else
                fhaplo.foundergeno[1][1,2] = ["1","N"]
                fhaplo.foundergeno[1][2,3] = ["N","2"]
                fhaplo.foundergeno[1][3,4] = ["N","N"]
            end
            @test isfile(savegenodata("founderhaplo2.vcf",fhaplo))
            @test isa(MagicBase.formfhaplo("founderhaplo2.vcf",
                isfounderinbred=isinbred), MagicGeno)
            @test isfile(savegenodata("founderhaplo2.csv",fhaplo))
            @test isa(MagicBase.formfhaplo("founderhaplo2.csv",
                isfounderinbred=isinbred), MagicGeno)
        end
        rm("founderhaplo2.vcf")
        rm("founderhaplo2.csv")
    end    
    verbose &&  println("4")
    @testset "magicsimulate_output" begin
        for dataid = ["simarray_magicsimulate","simmix_magicsimulate","multipop_magicsimulate"]
            filels  = [dataid*"_"*i*".csv.gz" for i=["truecontfgl",
                "truefgl","truegeno"]]
            pedfile = dataid*"_ped.csv"
            for fn= filels
                # @test isfile(fn)
                @test isa(formmagicgeno(fn,pedfile), MagicGeno)
                magicgeno = formmagicgeno(fn,pedfile)
                @test isa(savemagicgeno("magicgeno.csv",magicgeno), String)
                @test isa(readmagicgeno("magicgeno.csv"),MagicGeno)
            end
        end
        rm("magicgeno.csv")
    end    
    verbose &&  println("5")
    @testset "MagicAncestry" begin
        for dataid = ["simarray","simmix","multipop"]
            ancestryfile = dataid*"_magicreconstruct_ancestry.csv.gz"
            truefglfile= dataid*"_magicsimulate_truefgl.csv.gz"
            # @test isfile(ancestryfile)
            @test isa(readmagicancestry(ancestryfile),MagicAncestry)
            magicancestry = readmagicancestry(ancestryfile)
            @test isfile(savemagicancestry("magicancestry.csv",magicancestry))
            @test isa(readmagicancestry("magicancestry.csv"),MagicAncestry)
            pedfile = dataid*"_magicsimulate_ped.csv"            
            truegeno =formmagicgeno(dataid*"_magicsimulate_truegeno.csv.gz",pedfile)
            @test isa(formmagicgeno(dataid*"_magicsimulate_geno.vcf.gz",pedfile),MagicGeno)
            magicgeno = formmagicgeno(dataid*"_magicsimulate_geno.vcf.gz",pedfile)
            @test isa(MagicBase.rawgenoprob!(magicgeno), MagicGeno)
            @test isa(MagicBase.rawgenocall!(magicgeno), MagicGeno)            
        end
    end
    verbose && println("6")
    @testset "visualize" begin
        for dataid = ["simarray"]
            ancestryfile = dataid*"_magicreconstruct_ancestry.csv.gz"
            truefglfile= dataid*"_magicsimulate_truefgl.csv.gz"
            pedfile = dataid*"_magicsimulate_ped.csv"
            magicancestry = readmagicancestry(ancestryfile)
            truefgl = formmagicgeno(truefglfile,pedfile)
            @test isa(plotcondprob(magicancestry;truefgl,offspring=1),Any)
            @test isa(animcondprob(magicancestry;truefgl),Any)
        end
    end    
    rm("magicancestry.csv")
end




