

function reset_juncdist(juncdist::JuncDist,model::AbstractString;
    io::Union{IO,Nothing}=nothing,verbose::Bool=true, 
    isfounderinbred::Bool=true, isinferjunc::Bool=false)
    MagicBase.check_juncdist(juncdist)
    (; nfounder, ibd, mapexpansion, j1122, j1211, j1213, j1222,j1232) = juncdist
    b = isnothing.([j1122, j1211, j1213, j1222,j1232])
    if isnothing(mapexpansion)
        if any(b)
            if !isinferjunc                    
                msg = string("juncdensign is incomplete. juncdensign = ", juncdist)
                if model == "depmodel"
                    msg *= string(". JuncDist field: mapexpansion is required for model =",model)
                else
                    msg *= string(". JuncDist fields: ibd, j1122, j1211, j1213, j1222,j1232 are required for model =",model)
                end
                @error msg
                printconsole(io, false, msg)                       
            end
        else
            Rm = j1122+2*j1222+j1232
            Rp = j1122+2*j1211+j1213
            mapexpansion = (Rm + Rp)/2                                                
        end
    end    
    if model == "depmodel"
        ibd = 1.0
        isnothing(mapexpansion) && isinferjunc && (mapexpansion = 2.0)
        if !isnothing(mapexpansion)
            j1122 = mapexpansion
            j1211 = j1213 = j1222 = j1232 = 0.0
        end
    elseif model == "indepmodel"
        nfgl = isfounderinbred ? nfounder : 2*nfounder
        if isnothing(ibd) 
            ibd = 1.0/nfgl
        end
        isnothing(mapexpansion) && isinferjunc && (mapexpansion = 2.0)
        if !isnothing(mapexpansion)
            # j1232 = 0 and j1213 = 0 if nfgl = 2
            # Rm = j1122+2*j1222+j1232
            # Rp = j1122+2*j1211+j1213       
            R = mapexpansion             
            j1122 = ibd*R
            j1222 = j1211 = ((1-ibd)/nfgl)*R
            j1232 = j1213 = (1-ibd)*(1-2/nfgl)*R
            j1122, j1211, j1213, j1222, j1232 = round.((j1122, j1211, j1213, j1222,j1232),digits=4)
        end
    else
        model in ["jointmodel"] || @error string("unexpected model=",model) 
        if isnothing(ibd) 
            if isinferjunc                    
                ibd = 0.8
            else
                msg = string("juncdensign ibd is not provided. juncdensign = ", juncdist)
                @error msg
                printconsole(io, false, string("Error: ", msg))
            end
        else
            if any(b)
                if isnothing(mapexpansion)      
                    if isinferjunc           
                        mapexpansion = 2.0         
                    else
                        msg = string("juncdensign is incomplete. juncdensign = ", juncdist)
                        @error msg
                        printconsole(io, false, string("Error: ", msg))
                    end            
                end
                if !isnothing(mapexpansion)    
                    nfgl = isfounderinbred ? nfounder : 2*nfounder
                    R = mapexpansion             
                    j1122 = ibd*R
                    j1222 = j1211 = ((1-ibd)/nfgl)*R
                    j1232 = j1213 = (1-ibd)*(1-2/nfgl)*R
                    j1122, j1211, j1213, j1222, j1232 = round.((j1122, j1211, j1213, j1222,j1232),digits=4)
                end
            end
        end   
    end
    newjuncdist = JuncDist(; nfounder, ibd, mapexpansion,j1122, j1211, j1213, j1222,j1232)
    MagicBase.check_juncdist(newjuncdist)
    if isinferjunc           
        msg = string("set initial juncdist=", newjuncdist, " for ", model, 
            " and isinferjunc = true")
    else
        msg = string("reset juncdist=", newjuncdist, " for ", model)
    end
    printconsole(io, verbose, msg)   
    newjuncdist
end

function reset_juncdist!(magicped::MagicPed,model::AbstractString;
    io::Union{IO,Nothing}=nothing,verbose::Bool=true, 
    isfounderinbred::Bool=true, isinferjunc::Bool=false)
    if isa(magicped.designinfo, Dict{String,DesignInfo})
        for (popid,subdesign) in magicped.designinfo            
            if subdesign.designtype == :juncdist
                juncdist = parsedesign(subdesign.designcode).juncdist
                subdesign.juncdist = reset_juncdist(juncdist,model;io,verbose,isfounderinbred,isinferjunc) 
            end
        end        
    end
end

function reset_model(magicped::MagicPed,model::AbstractVector;
    io::Union{IO,Nothing},verbose::Bool)
    [reset_model(magicped,i;io, verbose) for i in model]
end

function reset_model(magicped::MagicPed,model::AbstractString;
    io::Union{IO,Nothing},verbose::Bool)
    offinfo = magicped.offspringinfo
    eltype(offinfo[!,:ishomozygous]) == Bool || @error string("magicped.offspringinfo[!,:ishomozygous] is not Bool")
    gdf = groupby(offinfo, :member)
    ishomozygousls = [begin
        ishomozygous = unique(df[!,:ishomozygous])
        if length(ishomozygous) > 1
            @error string("homozygous and heterozygous offspring are mixed in subpopulation=",df[1,:member])
        end
        only(ishomozygous)
    end for df in gdf]
    model2 = lowercase(model)
    if all(ishomozygousls)
        if model2 != "depmodel"
            model2 = "depmodel"
            msg = string("reset model=",model2," since all offspring are homozygous")
            MagicBase.printconsole(io,false,string("Warning: ", msg))
            verbose && @warn msg
        end
    end
    model2
end
