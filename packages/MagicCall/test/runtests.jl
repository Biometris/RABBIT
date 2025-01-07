using MagicCall
using MagicBase
using Test
cd(@__DIR__)

@testset "MagicCall" begin    
    @testset "outbred" begin
        @time include("2star-self1_outbred.jl")        
    end    
end
