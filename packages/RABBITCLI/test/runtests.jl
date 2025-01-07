using RABBITCLI
using Test
cd(@__DIR__)

@info string("RABBITCLI testing dir = ",pwd())

try 
    @static if Sys.isunix()        
        if !isfile("pipeline.sh")
            @warn string("RABBITCLI testing dir = ", pwd(), ", pipeline.sh does not exist!")
        end        
        run(`bash "pipeline.sh"`)
    elseif Sys.iswindows()
        if !isfile("pipeline.cmd")
            @warn string("RABBITCLI testing dir = ", pwd(), ", pipeline.cmd does not exist!")
        end                
        run(`"pipeline.cmd"`)        
    else
        @warn string("TODO: test RABBITCLI for ",Sys.MACHINE)
    end
    @test true
catch err
    @error string(err)        
    @test false
finally 
    # clean up
    cd(@__DIR__)
    rm.(filter(x->occursin("sim_", x), readdir()))
    rm.(filter(x->occursin(r"^jl_", x), readdir()))
end

