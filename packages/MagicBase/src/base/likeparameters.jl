
"""
    LikeParameters

keyword-based struct for the parameters of likelihood function. 

LikeParameters() is equivalent to LikeParameters(foundererror=0.005, offspringerror=nothing, peroffspringerror=0.0, 
seqerror=nothing, allelebalancemean=nothing, allelebalancedisperse=nothing, alleledropout=0.0). 
The `peroffspringerror` refers to error rate per offspring, and the other parameters refer to error rate per marker. 

If genotype format is not "AD", the parameters `seqerror`, `allelebalancemean`, `allelebalancedisperse`, and `alleledropout` are irrelevant. 

If `model="depmodel"`, the parameters `allelebalancemean`, `allelebalancedisperse`, and `alleledropout` are irrelevant. 

If there exists keyarg `isinfererror` and `isinfererror` = true, the parameters 
with values being nothing will be inferred and the other parameters will be fixed.  
If `isinfererror` = false,  LikeParameters() is equivalent to LikeParameters(foundererror=0.005, offspringerror=0.005, peroffspringerror=0.0, 
seqerror=0.001, allelebalancemean=0.5, allelebalancedisperse=0.0, alleledropout=0.0). 

"""
Base.@kwdef struct LikeParameters    
    foundererror::Union{Nothing,Float64} = 0.005    
    offspringerror::Union{Nothing,Float64} = nothing    
    peroffspringerror::Union{Nothing,Float64} = 0.0    
    seqerror::Union{Nothing,Float64} = nothing    
    allelebalancemean::Union{Nothing,Float64} = nothing    
    allelebalancedisperse::Union{Nothing,Float64} = nothing    
    alleledropout::Union{Nothing,Float64} = 0.0
end

function get_offspringerror(likeparameters::LikeParameters)
    (;offspringerror) = likeparameters
    isnothing(offspringerror) && (offspringerror = 0.005)
    offspringerror
end

function get_peroffspringerror(likeparameters::LikeParameters)
    (;peroffspringerror) = likeparameters
    isnothing(peroffspringerror) && (peroffspringerror = 0.0)
    peroffspringerror
end


function get_foundererror(likeparameters::LikeParameters)
    (;foundererror) = likeparameters
    isnothing(foundererror) && (foundererror = 0.005)
    foundererror
end

function get_seqerror(likeparameters::LikeParameters)
    (;seqerror) = likeparameters
    isnothing(seqerror) && (seqerror = 0.001)
    seqerror
end

function get_liketargetls(likeparameters::LikeParameters)    
    liketargetls  = String[]
    for i in propertynames(likeparameters)
        isnothing(getfield(likeparameters,i)) && push!(liketargetls,string(i))
    end
    liketargetls
end

function extract_likeparameters(likeparameters::LikeParameters)
    (; foundererror, offspringerror, peroffspringerror, seqerror, allelebalancemean,allelebalancedisperse,alleledropout) = likeparameters
    liketargetls  = String[]
    for i in propertynames(likeparameters)
        isnothing(getfield(likeparameters,i)) && push!(liketargetls,string(i))
    end		
    epsf = ifelse(isnothing(foundererror), 0.005, foundererror)
    epso = ifelse(isnothing(offspringerror), 0.005, offspringerror)
    epso_perind = ifelse(isnothing(peroffspringerror), 0.0, peroffspringerror)
    isnothing(seqerror) && (seqerror = 0.001)
    isnothing(allelebalancemean) && (allelebalancemean = 0.5)
    isnothing(allelebalancedisperse) && (allelebalancedisperse = 0.0)
    isnothing(alleledropout) && (alleledropout = 0.0)
    liketargetls, epsf, epso, epso_perind, seqerror, allelebalancemean,allelebalancedisperse,alleledropout
end

function parse_likeparameters(like::AbstractString)
    str = strip(like)
    if !occursin(r"^LikeParameters(.*)$",str)  
        msg = string("Could not parse likeparameters = ", str, "; it must be in format of LikeParameters()")
        error(msg)
    end    
    str = replace(str, "LikeParameters("=>"",")"=>"")
    isempty(str) && return LikeParameters()
    strls = strip.(split(str, ","))
    idset = ["foundererror", "offspringerror", "peroffspringerror","seqerror","allelebalancemean","allelebalancedisperse","alleledropout"]
    pairls = []    
    for i in strls
        idval = strip.(split(i,"="))
        length(idval) == 2 || error(string("Could not parse parameter: ",i))
        id, val = idval
        in(id,idset) || error(string("unknown parameter=",id,"; must be foundererror, offspringerror, peroffspringerror, seqerror, allelebalancemean,allelebalancedisperse, or alleledropout"))
        val2 = val == "nothing" ?  nothing : parse(Float64,val)
        push!(pairls,id=>val2)
    end    
    LikeParameters(; (Symbol(k) => v for (k,v) in pairls)...)    
end

"""
    ThreshLikeParameters

keyword-based struct for the thresholds of the likelihood parameters. Markers with the inferred parameter values 
being greater than the maximum will be deleted. 

ThreshLikeParameters() is equivalent to ThreshLikeParameters(foundererror=0.25, offspringerror=0.25, 
peroffspringerror=0.25, seqerror=0.25, allelebalancemean=0.9, allelebalancedisperse=1.0, alleledropout=0.05). 

"""
Base.@kwdef struct ThreshLikeParameters    
    foundererror::Float64 = 0.25
    offspringerror::Float64 = 0.25
    peroffspringerror::Float64 = 0.25
    seqerror::Float64 = 0.25
    allelebalancemean::Float64 = 0.9
    allelebalancedisperse::Float64 = 1.0
    alleledropout::Float64 = 0.05
end

function parse_threshlikeparameters(maxlike::AbstractString)
    str = strip(maxlike)
    if !occursin(r"^ThreshLikeParameters(.*)$",str)  
        msg = string("Could not parse threshlikeparameters = ", str, "; it must be in format of ThreshLikeParameters()")
        error(msg)
    end    
    str = replace(str, "ThreshLikeParameters("=>"",")"=>"")
    isempty(str) && return ThreshLikeParameters()
    strls = strip.(split(str, ","))
    idset = ["foundererror", "offspringerror", "peroffspringerror", "seqerror","allelebalancemean","allelebalancedisperse","alleledropout"]
    pairls = []    
    for i in strls
        idval = strip.(split(i,"="))
        length(idval) == 2 || error(string("Could not parse parameter: ",i))
        id, val = idval
        in(id,idset) || error(string("unknown parameter=",id,"; must be foundererror, offspringerror, peroffspringerror, seqerror, allelebalancemean,allelebalancedisperse, or alleledropout"))
        val2 = parse(Float64,val)
        push!(pairls,id=>val2)
    end    
    ThreshLikeParameters(; (Symbol(k) => v for (k,v) in pairls)...)    
end


"""
    PriorLikeParameters

keyword-based struct for the priors of the likelihood parameters. 

# Fields

`foundererror::Distribution = Beta(1.0,1.0)`: prior distribution of founder error rate

`offspringerror::Distribution = Beta(1.0,1.0)`: prior distribution of offspring error rate

`peroffspringerror::Distribution = Beta(1.0,1.0)`: prior distribution of error rate per offspring 

`seqerror::Distribution = Beta(1.0,1.0)`: prior distribution of sequence base error rate

`allelebalancemean::Distribution = Beta(1.01,1.01)` prior distribution of allele balance mean

`allelebalancedisperse::Distribution = Exponential(0.5)`: prior distribution of allele balance overdispersion

`alleledropout::Distribution = Beta(1,19)`: prior distribution of allele dropout rate

"""
Base.@kwdef struct PriorLikeParameters        
    foundererror::Distribution = Beta(1.0,1.0)
    offspringerror::Distribution = Beta(1.0,1.0)
    peroffspringerror::Distribution = Beta(1.0,1.0)
    seqerror::Distribution = Beta(1.0,1.0)
    allelebalancemean::Distribution = Beta(1.01,1.01)
    allelebalancedisperse::Distribution = Exponential(0.5)
    alleledropout::Distribution = Beta(1,19)
end

function parse_priorlikeparameters(priorlike::AbstractString)
    str = strip(priorlike)
    if !occursin(r"^PriorLikeParameters(.*)$",str)  
        msg = string("Could not parse priorlikeparameters = ", str, "; it must be in format of PriorLikeParameters()")
        error(msg)
    end        
    eval(Meta.parse(priorlike))
    # str = strip(replace(str, "PriorLikeParameters"=>""))
    # length(str) >=2 && (str = str[2:end-1])
    # isempty(str) && return PriorLikeParameters()
    # strls = strip.(split(str, "),")) # add ) to avoid , in distributions
    # idset = ["foundererror", "offspringerror", "peroffspringerror", "seqerror","allelebalancemean","allelebalancedisperse","alleledropout"]
    # pairls = []    
    # for i in strls
    #     idval = strip.(split(i,"="))
    #     length(idval) == 2 || error(string("Could not parse parameter: ",i))
    #     id, val = idval
    #     val *= ")" # add ) were removed by split ",)" in strls
    #     in(id,idset) || error(string("unknown parameter=",id,"; must be foundererror, offspringerror, peroffspringerror, seqerror, allelebalancemean,allelebalancedisperse, or alleledropout"))
    #     # val2 = parse(Float64,val)
    #     val2 = eval(Meta.parse(val))
    #     push!(pairls,id=>val2)
    # end    
    # PriorLikeParameters(; (Symbol(k) => v for (k,v) in pairls)...)    
end


