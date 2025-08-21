
"""
    LikeParam

keyword-based struct for the parameters of likelihood function. 

LikeParam() is equivalent to LikeParam(foundererror=0.005, offspringerror=nothing, peroffspringerror=0.0, 
baseerror=nothing, allelicbias=nothing, allelicoverdispersion=nothing, allelicdropout=0.0). 
The `peroffspringerror` refers to error rate per offspring, and the other parameters refer to error rate per marker. 

If genotype format is not "AD", the parameters `baseerror`, `allelicbias`, `allelicoverdispersion`, and `allelicdropout` are irrelevant. 

If `model="depmodel"`, the parameters `allelicbias`, `allelicoverdispersion`, and `allelicdropout` are irrelevant. 

If there exists keyarg `isinfererror` and `isinfererror` = true, the parameters 
with values being nothing will be inferred and the other parameters will be fixed.  
If `isinfererror` = false,  LikeParam() is equivalent to LikeParam(foundererror=0.005, offspringerror=0.005, peroffspringerror=0.0, 
baseerror=0.001, allelicbias=0.5, allelicoverdispersion=0.0, allelicdropout=0.0). 

"""
Base.@kwdef struct LikeParam    
    foundererror::Union{Nothing,Float64} = 0.005    
    offspringerror::Union{Nothing,Float64} = nothing    
    peroffspringerror::Union{Nothing,Float64} = 0.0    
    baseerror::Union{Nothing,Float64} = 0.001    
    allelicbias::Union{Nothing,Float64} = nothing    
    allelicoverdispersion::Union{Nothing,Float64} = nothing    
    allelicdropout::Union{Nothing,Float64} = 0.0
end

function get_likeproperty(likeparam::LikeParam, id)
    (;foundererror,offspringerror,peroffspringerror,baseerror,allelicbias,allelicoverdispersion,allelicdropout) = likeparam
    if id == :foundererror
        isnothing(foundererror) ? 0.005 : foundererror
    elseif id == :offspringerror
        isnothing(offspringerror) ? 0.005 : offspringerror 
    elseif id == :peroffspringerror
        isnothing(peroffspringerror) ? 0.0 : peroffspringerror
    elseif id == :baseerror
        isnothing(baseerror) ? 0.001 : baseerror
    elseif id == :allelicbias
        isnothing(allelicbias) ? 0.5 : allelicbias
    elseif id == :allelicoverdispersion
        isnothing(allelicoverdispersion) ? 0.0 : allelicoverdispersion
    elseif id == :allelicdropout
        isnothing(allelicdropout) ? 0.0 : allelicdropout
    else
        @error string("unknown id=",id,"; must be foundererror, offspringerror, peroffspringerror, baseerror, allelicbias, allelicoverdispersion, or allelicdropout")        
    end
end

function get_liketargetls(likeparam::LikeParam)    
    liketargetls  = String[]
    for i in propertynames(likeparam)
        isnothing(getfield(likeparam,i)) && push!(liketargetls,string(i))
    end
    liketargetls
end

function extract_likeparam(likeparam::LikeParam)
    (; foundererror, offspringerror, peroffspringerror, baseerror, allelicbias,allelicoverdispersion,allelicdropout) = likeparam
    liketargetls  = String[]
    for i in propertynames(likeparam)
        isnothing(getfield(likeparam,i)) && push!(liketargetls,string(i))
    end		
    epsf = ifelse(isnothing(foundererror), 0.005, foundererror)
    epso = ifelse(isnothing(offspringerror), 0.005, offspringerror)
    epso_perind = ifelse(isnothing(peroffspringerror), 0.0, peroffspringerror)
    isnothing(baseerror) && (baseerror = 0.001)
    isnothing(allelicbias) && (allelicbias = 0.5)
    isnothing(allelicoverdispersion) && (allelicoverdispersion = 0.0)
    isnothing(allelicdropout) && (allelicdropout = 0.0)
    liketargetls, epsf, epso, epso_perind, baseerror, allelicbias,allelicoverdispersion,allelicdropout
end


"""
    SoftThreshLikeParam

keyword-based struct for the thresholds of the likelihood parameters. 
In magiccall and magicimpute, markers with the inferred parameter values being greater than the thresholds will be deleted only if they are also outerliers. 
In magicimpue, offspring with the inferred parameter values being greater than the thresholds will be excluded only if they are also outerliers. 

SoftThreshLikeParam() is equivalent to SoftThreshLikeParam(foundererror=0.025, offspringerror=0.025, 
peroffspringerror=0.025, baseerror=0.01, allelicbias=0.67, allelicoverdispersion=0.5, allelicdropout=0.01). 

"""
Base.@kwdef struct SoftThreshLikeParam    
    foundererror::Float64 = 0.025
    offspringerror::Float64 = 0.025
    peroffspringerror::Float64 = 0.025
    baseerror::Float64 = 0.01
    allelicbias::Float64 = 0.67
    allelicoverdispersion::Float64 = 0.5
    allelicdropout::Float64 = 0.01
end

"""
    ThreshLikeParam

keyword-based struct for the thresholds of the likelihood parameters. 
In magiccall and magicimpute, markers with the inferred parameter values being greater than the thresholds will be deleted. 
In magicimpue, offspring with the inferred parameter values being greater than the thresholds will be excluded. 

ThreshLikeParam() is equivalent to ThreshLikeParam(foundererror=0.25, offspringerror=0.25, 
peroffspringerror=0.25, baseerror=0.05, allelicbias=0.9, allelicoverdispersion=1.0, allelicdropout=0.1). 

"""
Base.@kwdef struct ThreshLikeParam    
    foundererror::Float64 = 0.25
    offspringerror::Float64 = 0.25
    peroffspringerror::Float64 = 0.25
    baseerror::Float64 = 0.05
    allelicbias::Float64 = 0.9
    allelicoverdispersion::Float64 = 1.0
    allelicdropout::Float64 = 0.1
end

"""
    PriorLikeParam

keyword-based struct for the priors of the likelihood parameters. 

# Fields

`foundererror::Union{Nothing,Distribution} = nothing`: prior distribution of founder error rate. 
If nothing, it is given by Beta(1,1/merr-1) where merr is the mean among markers in the previous iteration.   

`offspringerror::Union{Nothing,Distribution} = nothing`: prior distribution of offspring error rate
If nothing, it is given by Beta(1,1/merr-1) where merr is the mean among markers in the previous iteration.   

`peroffspringerror::Union{Nothing,Distribution} = nothing`: prior distribution of error rate per offspring 
If nothing, it is given by Beta(1,1/merr-1) where merr is the mean among offspring in the previous iteration.   

`baseerror::Union{Nothing,Distribution} = nothing`: prior distribution of sequence base error rate
If nothing, it is given by Beta(1,1/merr-1) where merr is the mean among markers in the previous iteration.   

`allelicbias::Union{Nothing,Distribution} = nothing`: prior distribution of allele balance mean. 
If nothing, it is given by Beta(1.01,1.01). 

`allelicoverdispersion::Union{Nothing,Distribution} = nothing`: prior distribution of allele balance overdispersion
If nothing, it is given by Exponential(merr) where merr is the mean among markers in the previous iteration.   

`allelicdropout::Union{Nothing,Distribution} = nothing`: prior distribution of allele dropout rate
If nothing, it is given by Beta(1,1/merr-1) where merr is the mean among markers in the previous iteration.   

"""
Base.@kwdef struct PriorLikeParam        
    foundererror::Union{Nothing,Distribution} = nothing
    offspringerror::Union{Nothing,Distribution} = nothing
    peroffspringerror::Union{Nothing,Distribution} = nothing
    baseerror::Union{Nothing,Distribution} = nothing
    allelicbias::Union{Nothing,Distribution} = nothing
    allelicoverdispersion::Union{Nothing,Distribution} = nothing
    allelicdropout::Union{Nothing,Distribution} = nothing
end

function get_info_likeparam(likeparam::Union{LikeParam,SoftThreshLikeParam, ThreshLikeParam,PriorLikeParam}; 
    isinfererror=true, ismultiline=false)
    msgls = []
    t = typeof(likeparam)
    push!(msgls,string(t,"("))
    namels = propertynames(likeparam)
    for id in namels
        val = getproperty(likeparam,id)
        if isnothing(val) && isa(likeparam,LikeParam)
            if isinfererror 
                val = "to_be_inferred" 
            else 
                val = get_likeproperty(likeparam, id) 
            end
        end
        push!(msgls,string(id," = ",val,","))
    end
    if ismultiline
        msg = join(msgls,"\n\t")
        msg = string(rstrip(msg,','), "\n\t)")
    else
        msg = join(msgls,"")
        msg = string(rstrip(msg,','), ")")
    end    
    msg
end

function parse_likeparam(like::AbstractString, type::AbstractString)
    str = strip(like)
    str == "nothing" && return nothing    
    if !occursin(Regex(string("^",type,"(.*)\$")),str)  
        msg = string("Could not parse ", str, ": it must be in format of "*type*"()")
        error(msg)
    end        
    eval(Meta.parse(like))
end
