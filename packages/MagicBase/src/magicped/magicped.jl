
"""
    MagicPed

mutable struct that stores pedigree information

    MagicPed(designinfo,founderinfo,offspringinfo)

inner constructor. See also [`readmagicped`](@ref).

# Fields

`designinfo::Union{Nothing, Dict{String,DesignInfo},Pedigree}`: specifies population
designinfo. A designinfo::Pedigree specifies the designinfo in pedigree. See also [`Pedigree`](@ref). 
A designinfo::Dict{String,DesignInfo} specifies designinfo for each subpopulation. 
See also [`DesignInfo`](@ref). 

`founderinfo::DataFrame`: founder information. 

`offspringinfo::DataFrame`: offspring information. The column names are
[:individual,:member,:ishomozygous,:isfglexch,:gender]. The :ishomozygous column specifies if
the individual is homozygous. If design is set by juncdist::JuncDist,
the :member column is set to "juncdist". If design is set by a string
designcode::AbstractString, the :member column is set to the last non-male
member of the pedigree. If design is set by designinfo::Pedigree,
the :member column is associated to the :member column of the pedigree.

"""
mutable struct MagicPed
    designinfo::Union{Nothing,Dict{String,DesignInfo}, Pedigrees.Pedigree}
    #founderinfo [:individual,:gender]
    founderinfo::Union{Nothing,DataFrame}
    #offspringinfo [:individual,:member, :ishomozygous,:isfglexch,:gender]
    offspringinfo::Union{Nothing,DataFrame}
    function MagicPed(designinfo,founderinfo,offspringinfo)
        new(designinfo,founderinfo,offspringinfo)
    end
end


function formmagicped(genofile::AbstractString,pedinfo::AbstractString;    
    ishomozygous::Bool=false, 
    delim::AbstractChar=',',
    commentstring::AbstractString="##",    
    workdir::AbstractString=pwd())
    if isa(pedinfo, AbstractString)
        if last(splitext(pedinfo))==".csv"
            # pedfile
            magicped = readmagicped(pedinfo; delim, commentstring, workdir)        
        else
            # designcode
            designinfo = MagicBase.parsedesign(pedinfo)                 
            inext = last(MagicBase.split_allext(genofile))
            if inext in [".vcf",".vcf.gz"] 
                samplels = vcf_get_samples(genofile; commentstring, workdir)
            elseif inext in [".csv",".csv.gz"] 
                samplels = csv_get_samples(genofile; commentstring, workdir)
            else
                @error string("genofile=",genofile, " has unknown ext = ", inext)
            end            
            magicped = formmagicped!(designinfo,samplels; ishomozygous)    
        end
    else
        # juncdist        
        @error string("TODO: formmagicped for pedinfo=",pedinfo)
    end
end


"""
    formmagicped(designinfo, founderids,offspringids)

form magicped::MagicPed for a non-divided pouplation with indiviudal ids. 

"""
function formmagicped!(designinfo::DesignInfo,
    samplels::AbstractVector;
    isfglexch::Union{Nothing,Bool}=nothing,
    ishomozygous::Bool=false)    
    if designinfo.designtype == :breedcross
        founderids = designinfo.founders
        issubset(founderids,samplels) || error(string("designinfo.founders = ",founderids," are not in samplels"))
    else
        nf = getnfounder(designinfo)
        founderids = samplels[1:nf]
    end
    offspringids = setdiff(samplels,founderids)
    popsize = length(offspringids)
    magicped = formmagicped(designinfo,popsize; ishomozygous,isfglexch) 
    setfounderid!(magicped,founderids)
    setoffspringid!(magicped,offspringids)
    magicped
end


"""
    formmagicped(designinfo, popsize)

form magicped::MagicPed form magicped for a non-divided pouplation of popsize. 

"""
function formmagicped(designinfo::DesignInfo,popsize::Integer;
    isfglexch::Union{Nothing,Bool}=nothing,
    ishomozygous::Bool=false)
    if designinfo.designtype == :commoncross
        if isnothing(isfglexch)
            isfglexch = occursin(r"^[1-9][0-9]{0,6}ril-",designinfo.designcode)                
        end
        formmagicped(designinfo.pedigree,popsize;ishomozygous,isfglexch)
    elseif designinfo.designtype == :breedcross
        if isnothing(isfglexch)
            isfglexch = false
        end
        formmagicped(designinfo.pedigree,popsize;ishomozygous,isfglexch)
    elseif designinfo.designtype == :juncdist
        formmagicped(designinfo.juncdist,popsize; popid="subpop") # isfglexch=true
    else
        error(string("unknown designtype=",designinfo.designtype))
        0
    end 
end

"""
    formmagicped(pedigree, popsize)

form magicped::MagicPed from pedigree and popsize. 

"""
function formmagicped(pedigree::Pedigrees.Pedigree,popsize::Integer;
    ishomozygous::Bool=false,isfglexch::Bool=false)
    # offspringinfo 
    subpopls = findall(pedigree.generation .== last(pedigree.generation))
    # isnothing(popsize) && (popsize = length(subpopls))
    popsize>0 || @error string("popsize=",popsize, " is not positive")
    if length(subpopls) > popsize
        ii = sort(sample(1:length(subpopls),popsize,replace=false))
        subpopls = subpopls[ii]
    end
    nsubpop = length(subpopls)
    persize,rem = divrem(popsize, nsubpop)
    sizels = persize*ones(Int, nsubpop)
    rem>0 && (sizels[1:rem] .+= 1)    
    offspringinfo =vcat([DataFrame(:individual=>[string("offspring_",i,"_",j) for j=1:sizels[i]],
        :member=>pedigree.member[subpopls[i]],
        :ishomozygous=>ishomozygous,
        :isfglexch=>isfglexch,
        :gender=>pedigree.gender[subpopls[i]]) for i=1:nsubpop]...)
    # founderinfo 
    nf = pedigree.nfounder
    founderinfo = DataFrame(:individual=>pedigree.member[1:nf],
        :gender=>pedigree.gender[1:nf])
    MagicPed(pedigree,founderinfo,offspringinfo)
end

"""
    formmagicped(designcode, popsize)

form magicped::MagicPed from designcode and popsize. 

"""
function formmagicped(designcode::AbstractString,popsize::Integer)
    designinfo = parsedesign(designcode)
    formmagicped(designinfo,popsize)
end

"""
    formmagicped(juncdist, popsize)

form magicped::MagicPed from juncdist and popsize. 

"""
function formmagicped(juncdist::JuncDist,popsize::Integer; popid::AbstractString="subpop")
    popsize>=0 || @error string("popsize=",popsize, " is negative")
    check_juncdist(juncdist)
    #offspringinfo
    ishomozygous = juncdist.ibd == 1.0
    offspringinfo =DataFrame(:individual=>[string("offspring_",j) for j=1:popsize],
        :member=>popid,
        :ishomozygous=>ishomozygous,
        :isfglexch => true,
        :gender=>"notapplicable")
    # founderinfo 
    nf = juncdist.nfounder
    founderinfo = DataFrame(:individual=>[string("founder_",i) for i in 1:nf],
        :gender=>"notapplicable")
    designcode = todesigncode(juncdist)
    designinfo = Dict(popid=>DesignInfo(; designtype=:juncdist, designcode, juncdist))
    MagicPed(designinfo,founderinfo,offspringinfo)
end

function setfounderid!(magicped::MagicPed,founders::AbstractVector)    
    nf = size(magicped.founderinfo,1)
    if length(founders) < nf
        @error string("the number of input founders must be  ",nf)
    elseif length(founders) > nf
        founders = founders[1:nf]
    end
    founders=string.(founders)  
    if isa(magicped.designinfo,Pedigrees.Pedigree)
        magicped.designinfo = Pedigrees.setfounderid(magicped.designinfo,founders)
    elseif isa(magicped.designinfo,Dict{String,DesignInfo})
        if length(magicped.designinfo) == 1
            design = last(first(magicped.designinfo))
            setfounders!(design,founders)     
        else
            @error string("Could not setfounderid for designinfo dict length =",length(magicped.designinfo))
        end
    end
    magicped.founderinfo[!,:individual] = founders
    magicped
end


function setoffspringid!(magicped::MagicPed,offspring::AbstractVector)
    noff = size(magicped.offspringinfo,1)
    if length(offspring) < noff
        error("the number of input offspring must be ",noff)
    elseif length(offspring) > noff
        offspring = offspring[1:nf]
    end
    magicped.offspringinfo[!,:individual] = string.(offspring)
    magicped
end


function df2designinfo(designdf::DataFrame)
    isempty(designdf) && return nothing
    designcols = propertynames(designdf)
    in(:member, designcols) || @error string("pedfile designinfo must have :member column")
    allunique(designdf[!,:member]) || @error string("pedfile designinfo :member column must be allunique")    
    if issubset([:mother,:father],designcols)
        designinfo = Pedigrees.Pedigree(designdf)
    elseif in(:designcode,designcols)
        hasfounders = in(:founders,designcols)
        designinfo = Dict([begin 
            if hasfounders
                if ismissing(row[:founders])
                    founders0 = "NA"
                    founders = nothing
                else
                    founders0 = strip(row[:founders])
                    for sep in ["/",",","=","=>"]
                        if occursin(sep, founders0) 
                            @error string(sep, " is not expected in founders=", founders0, "; sperate founders by ||")
                        end
                    end
                    founders = string.(strip.(split(founders0,"||")))
                end
            else
                founders0 = "NA"
                founders = nothing
            end
            popid = convert(String,row[:member])          
            designcode = row[:designcode]    
            design = parsedesign(designcode; founders,popid)    
            if design.designtype == :breedcross 
                if !in(founders0, ["NA","missing"]) 
                    if !issetequal(strip.(split(founders0,"||")), design.founders)
                        @warn string("founders=", founders0, ", inconsistent with parsed founders=",design.founders)
                    end
                end
            else
                in(founders0, ["NA","missing"]) && @warn string("founders=", founders0, " are missing for designcode=", esigncode, " with designtype=",design.designtype)
            end
            popid => design
        end for (i,row) in enumerate(eachrow(designdf))])
    else
        @error string("Could not parse pedfile designinfo columns: ",designcols, "; either has columns [:mother, :father] or has columns [:designcode,:founders]")
    end
    designinfo
end

function parsemagicped(peddict::AbstractDict)
    strkeys = collect(keys(peddict))
    if !issubset(["designinfo","offspringinfo"],strkeys)
        error(string("keys designinfo and offspringinfo are not in ",strkeys))
    end    
    # desigininfo    
    designdf = peddict["designinfo"]
    designinfo = df2designinfo(designdf)   
    # offspringinfo
    offspringinfo = peddict["offspringinfo"]
    if !isempty(offspringinfo)        
        checkoffspringinfo!(offspringinfo,designinfo)
    end
    founderinfo = get_founderinfo(designinfo)    
    MagicPed(designinfo,founderinfo, offspringinfo)    
end

function checkoffspringinfo!(offspringinfo::DataFrame,designinfo::Union{Nothing,Dict{String,DesignInfo}, Pedigrees.Pedigree})
    offcols = propertynames(offspringinfo)
    in(:individual, offcols) || @error string("pedfile offspringinfo must have :individual column")
    in(:member, offcols) || @error string("pedfile offspringinfo must have :member column")
    for i in [:individual,:member]
        offspringinfo[!,i]=string.(strip.(string.(offspringinfo[!,i])))
    end
    memls = unique(offspringinfo[!,:member])            
    if isa(designinfo,Pedigrees.Pedigree)
        d = setdiff(memls,designinfo.member)
        isempty(d) || @error string("members in offspringinfo but not in designinfo pedigree: ",d)
    elseif isa(designinfo, Dict{String, DesignInfo})
        d = setdiff(memls,keys(designinfo))
        isempty(d) || @error string("members in offspringinfo but not in designinfo dict: ",d)            
        d = setdiff(keys(designinfo),memls)
        isempty(d) || @warn string("members in designinfo dict but not in doffspringinfo: ",d)        
    else
        isnothing(designinfo) || @error string("unexpected designinfo=",designinfo)    
    end
    if !in(:ishomozygous, offcols)
        insertcols!(offspringinfo,3,:ishomozygous => falses(size(offspringinfo,1)))            
    else
        if eltype(offspringinfo[!,:ishomozygous]) <: AbstractString
            offspringinfo[!,:ishomozygous] .= [parse(Bool, lowercase(i)) for i in offspringinfo[!,:ishomozygous]]
        end
    end
    if !in(:isfglexch, offcols)
        insertcols!(offspringinfo,4,:isfglexch => falses(size(offspringinfo,1)))
    else
        if eltype(offspringinfo[!,:isfglexch]) <: AbstractString
            offspringinfo[!,:isfglexch] .= [parse(Bool, lowercase(i)) for i in offspringinfo[!,:isfglexch]]
        end
    end 
    for col in offcols
        if occursin(r"^peroffspringerror",string(col))
            if eltype(offspringinfo[!,col]) <: AbstractString
                offspringinfo[!,col] .= [parse(Float64, i) for i in offspringinfo[!,col]]
            end
        end
        if occursin(r"^offspringexcl",string(col))
            if eltype(offspringinfo[!,col]) <: AbstractString
                offspringinfo[!,col] .= [parse(Bool, i) for i in offspringinfo[!,col]]
            end
        end
    end
    if in(:isoutlier, offcols)       
        if eltype(offspringinfo[!,:isoutlier]) <: AbstractString
            offspringinfo[!,:isoutlier] .= [parse(Bool, lowercase(i)) for i in offspringinfo[!,:isoutlier]]
        end
    end 
    if !in(:gender, offcols)        
        if isa(designinfo, Pedigrees.Pedigree)
            dict = Dict(designinfo.member .=> designinfo.gender)
            genderls = [get(dict, i, nothing) for i in memls]
            dict = Dict(memls .=> genderls)
            offgender = [get(dict,i,missing) for i=offspringinfo[!,:member]]
            insertcols!(offspringinfo,5,:gender => offgender)
        else
            insertcols!(offspringinfo,5,:gender => "notapplicable")
        end
    end        
    memls = unique(offspringinfo[!,:member])            
    for mem in memls
        b = offspringinfo[!,:member] .== mem
        df = unique(offspringinfo[b,[:ishomozygous,:isfglexch,:gender]])
        if size(df,1) > 1            
            for col in [:ishomozygous,:isfglexch,:gender]
                length(unique(df[!,col])) > 1 || @error string("pop=",mem, " has non-unique ", col, ": ", unique(df[!,col]))
            end
        else
            if isa(designinfo, Dict{String, DesignInfo})                    
                designtype = designinfo[mem].designtype                    
                if designtype == :juncdist
                    isfglexch = df[1,:isfglexch] 
                    if !isfglexch
                        @info string("set isfglexch=true for pop=", mem," with designtype=",designtype)                            
                        offspringinfo[b,:isfglexch] .= true
                    end
                    ibd = designinfo[mem].juncdist.ibd
                    ishomo = df[1,:ishomozygous]
                    if ibd ≈ 1.0 && !ishomo
                        @info string("set ishomozygous=true for pop=",mem, " with ibd ≈ 1.0")                            
                        offspringinfo[b,:ishomozygous] .= true      
                    elseif ibd ≈ 0.0 && ishomo     
                        @warn string("set ishomozygous=false for pop=",mem, " with ibd ≈ 0.0")                            
                        offspringinfo[b,:ishomozygous] .= false                 
                    end
                end   
                designcode = designinfo[mem].designcode
                if occursin(r"=>DH$",designcode) && !df[1,:ishomozygous]
                    @info string("set ishomozygous=true for DH_pop=",mem)                            
                    offspringinfo[b,:ishomozygous] .= true    
                elseif occursin(r"=>FIXED",designcode) && !df[1,:ishomozygous]
                    @info string("set ishomozygous=true for FIXED_pop=",mem)                            
                    offspringinfo[b,:ishomozygous] .= true    
                end
            end
        end
    end   
    offspringinfo
end

function get_founderinfo(designinfo)
    if isa(designinfo,Pedigree)
        nf=designinfo.nfounder
        founderinfo = DataFrame([designinfo.member[1:nf] designinfo.gender[1:nf]],[:individual,:gender])
    elseif isa(designinfo,Dict{String,DesignInfo})
        founders = sort(unique(reduce(vcat, [i.founders for i in values(designinfo)])))
        founderinfo = DataFrame(individual=founders,gender="notapplicable")        
    else
        founderinfo = nothing
    end
    founderinfo
end

"""
    readmagicped(pedfile; commentstring="##",workdir=pwd())

read a CSV formatted `pedfile` and return magicped::MagicPed. 

# Positional arguments

`pedfile::AbstractString`: saves the pedigre information: designinfo and offspringinfo. 
The designinfo can be provided in three formats: pedigree, mating design code, 
and junction distribution. 
  
# Keyword arguments

`commentstring::AbstractString="##"`: the lines beginning with commentstring are ignored
in pedfile.

`workdir::AbstractString=pwd()`: directory for reading pedfile.

"""
function readmagicped(pedfile::AbstractString;
        delim::AbstractChar=',',
        commentstring::AbstractString="##",
        workdir::AbstractString=pwd())
    pedfile2=getabsfile(workdir,pedfile)    
    peddict = readmultitable(pedfile2; delim,commentstring,isparse=false)
    parsemagicped(peddict)
end

"""
    savemagicped(sink,magicped; delim=',',workdir=pwd())

save magicped into sink. see [`readmagicped`] for reading saved output file.

# Positional arguments

`sink::Union{IO,AbstractString}`: output file or IO.

`magicped::MagicPed`: a struct for storing pedidgree info.

# Keyword arguments

`delim::AbstractChar=','`: delimitor character.

`workdir::AbstractString=pwd()`: directory for reading genfile and pedfile.

"""
function savemagicped(sink::Union{IO,AbstractString},magicped::MagicPed;
    workdir::AbstractString = pwd(),
    initial::AbstractString="RABBIT", 
    delim::AbstractChar=',')
    if typeof(sink) <: AbstractString
        outputfile2 = getabsfile(workdir,sink)
        io = open(outputfile2, "w")
    elseif typeof(sink) <: IO
        io = sink
    else
        error(string("unknow sink: ",sink))
    end    
    ped=magicped.designinfo
    df = designinfo2df(ped)    
    appenddf(io, df; delim,initial,dfname="designinfo")
    # df cols: indiviudal, member, ishomozygous, gender    
    df = offspringinfo2df(magicped.offspringinfo)
    appenddf(io, df; delim,initial,dfname="offspringinfo")
    flush(io)
    typeof(sink) <: AbstractString && close(io)
    sink
end

function offspringinfo2df(offspringinfo)
    isnothing(offspringinfo) && return nothing
    df  = offspringinfo
    if !isnothing(df)
        cols = propertynames(df)
        if in(:gender,cols)
            isnagender = all(unique(df[!,:gender]) .== "notapplicable")        
            if isnagender
                setdiff!(cols, [:gender])
            end
        end        
        if in(:ishomozygous,cols) && !any(df[!,:ishomozygous])
            setdiff!(cols, [:ishomozygous])
        end
        if in(:isfglexch,cols) && !any(df[!,:isfglexch])
            setdiff!(cols, [:isfglexch])
        end
        df = df[!,cols]                   
    end
    df
end

function designinfo2df(designinfo::Union{Nothing,Dict{String,DesignInfo}, Pedigrees.Pedigree})
    if isnothing(designinfo)
        df = nothing
    else
        if isa(designinfo, Pedigree)
            title=[fieldnames(typeof(designinfo))[2:end]...]
            df=DataFrame([getfield(designinfo,i) for i=title],title)
            # df cols: member, mother, father, gender, generation
            isnagender = all(unique(df[!,:gender]) .== "notapplicable")
            isnagender && (df = df[!, 1:3])
            df
        elseif isa(designinfo,Dict{String,DesignInfo})
            # ids = reduce(vcat,propertynames(subdesign))
            # setdiff!(ids,[:nfounder,:founders])
            # juncvals = [getfield(juncdist,i) for i in ids]
            # b = .!isnothing.(juncvals)
            # juncstr = join(map((x,y)->string(x,"=",y),juncids[b],juncvals[b]),"||")    
            isbreedcross = [sub.designtype == :breedcross for (_,sub) in  designinfo]
            if all(isbreedcross)
                data = [[popid subdesign.designcode]    
                    for (popid,subdesign) in designinfo]
                df = DataFrame(reduce(vcat,data),[:member,:designcode])
            else
                data = [begin 
                    designcode = isnothing(subdesign.designcode) ? "NA" : subdesign.designcode
                    founders = subdesign.designtype == :breedcross ? "NA" : join(subdesign.founders,"||")
                    [popid designcode founders]    
                end for (popid,subdesign) in designinfo]
                df = DataFrame(reduce(vcat,data),[:member,:designcode,:founders])                
            end
            sort!(df,:member)
        else
            error(string("unknown designinfor type: ",typeof(designinfo)))
        end        
    end
end


"""
    plotmagicped(magicped; kwargs...)

plot magicped. 

# Keyword arguments

`isfounderinbred::Bool=true`: if true, founders are inbred, and otherwise outbred.

`outfile::Union{Nothing,AbstractString}=nothing`: if nothing, not save the plot, and otherwise save to outfile. 

"""
function plotmagicped(magicped::MagicPed;
    isfounderinbred::Bool=true,
    curves = false,
    names = nothing,
    nodecolor = nothing,
    plotsize = nothing,
    markersize = nothing,
    fontsize = 10,
    nodesize = 0.1,
    edgecolor = :black,
    edgestyle = :solid,    
    outfile::Union{Nothing,AbstractString}=nothing,
    plotkeyargs... 
    )
    infotype = typeof(magicped.designinfo) 
    if infotype <: Pedigree || infotype <: Dict{String,DesignInfo}
        ped =  magicped2ped(magicped;isfounderinbred)
        isnothing(ped) && return string("magicped.designinfo: ", magicped.designinfo)
        nf=ped.nfounder
        if isnothing(nodecolor)
            memls = magicped.offspringinfo[!,:member]
            offmem=unique(memls)
            bool=[i in offmem for i in ped.member]
            if isfounderinbred
                founder_colors = distinguishable_colors(nf, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
            else
                founder_colors = [RGB(0.9,0.9,0.9) for _ in 1:nf]
            end
            nodecolor = vcat(founder_colors,[RGB(0.9,0.9,0.9) for i=nf+1:length(ped.member)])
            nodecolor[bool] .= RGB(0.5,0.5,0.5)  
        end
        if isnothing(names)            
            memls = magicped.offspringinfo[!,:member]            
            offmem=unique(memls)
            nodename = [" " for _ in 1:length(ped.member)]
            nodename[1:nf] .=  ped.member[1:nf]
            for i in nf+1:length(ped.member)
                mem = ped.member[i]
                if in(mem,offmem)
                    b = memls .== mem
                    ishomozygousls = magicped.offspringinfo[b,:ishomozygous]
                    poptype = all(ishomozygousls) ? "Homo" : "N"
                    nodename[i] = string(mem,"\n",poptype,"=",sum(b))
                    if all(ishomozygousls) 
                        nodecolor[i] = RGB(0.7,0.7,0.7)
                        nodename[i] = string(mem,"\nHomo=",sum(b))
                    else
                        nodename[i] = string(mem,"\nN=",sum(b))
                    end
                end
            end
        else
            nodename = names
        end
        fig = Pedigrees.plotped(ped;names=nodename,nodecolor,
            markersize, plotsize,
            curves,fontsize,nodesize,edgecolor,edgestyle,
            plotkeyargs...
        )
        if !isnothing(outfile) 
            try  
                savefig(fig,outfile)
            catch
                @warn string("Cannot savefig to ",outfile)
            end
        end
        fig
    else
        string("magicped.designinfo: ", magicped.designinfo)
    end
end

function submagicped(magicped::MagicPed,subpop::AbstractString)
    submagicped(magicped,[subpop])
end

function submagicped(magicped::MagicPed,subpops::AbstractVector)
    memls = unique(magicped.offspringinfo[!,:member])
    issubset(subpops,memls) || error(string("subpops=",subpops, " is not a subset of subpopulations = ",memls))
    newoffinfo = filter(row->in(row[:member],subpops),magicped.offspringinfo)
    if isa(magicped.designinfo, Dict{String, DesignInfo})        
        newdesigninfo = Dict([convert(String,subpop)=>magicped.designinfo[subpop] for subpop in subpops])         
        founders = unique(reduce(vcat, [newdesigninfo[subpop].founders for subpop in subpops]))        
        newfounderinfo = filter(row->in(row[:individual],founders),magicped.founderinfo)        
        MagicPed(newdesigninfo,newfounderinfo,newoffinfo)
    elseif isa(magicped.designinfo, Pedigree)
        # designinfo        
        newdesigninfo = Pedigrees.getsubped(magicped.designinfo, subpops)
        # founderinfo        
        founders = Pedigrees.getfounderid(newdesigninfo)
        newfounderinfo = filter(row->in(row[:individual],founders),magicped.founderinfo)
        MagicPed(newdesigninfo,newfounderinfo,newoffinfo)
    else
        @error string("Could not extract ped from designinfo=",magicped.designinfo)
        return nothing        
    end
end


"""
    pedfile_designcode2ped(pedfile; commentstring='#',workdir=pwd())

convert a pedfile from `designinfo` being `designcode` to `pedigree`.

# Keyword arguments

`isfounderinbred::Bool=true`: if true, founders are inbred, and otherwise outbred. 

`commentstring::AbstractString="##"`: the lines beginning with commentstring are ignored
in pedfile.

`workdir::AbstractString=pwd()`: directory for reading pedfile.

"""
function pedfile_designcode2ped(pedfile::AbstractString;
    isfounderinbred::Bool=true,    
    outfile = "oustem.csv", 
    commentstring::AbstractString="##",
    workdir::AbstractString=pwd())
    magicped = MagicBase.readmagicped(pedfile;commentstring,workdir)
    magicped_designcode2ped!(magicped;isfounderinbred)
    MagicBase.savemagicped(outfile,magicped; workdir)
end

function magicped_designcode2ped!(magicped::MagicPed;
  isfounderinbred::Bool=true)
  ped = MagicBase.magicped2ped(magicped;isfounderinbred)  
  allfounders = ped.member[1:ped.nfounder]
  if magicped.founderinfo[!,:individual] != allfounders
    dict = Dict(magicped.founderinfo[!,:individual] .=> 1:size(magicped.founderinfo,1))
    ii = [dict[i] for i in allfounders]
    magicped.founderinfo = magicped.founderinfo[ii,:]
  end  
  popidls = []
  for (popid,subdesign) in magicped.designinfo
    in(subdesign.designtype,[:commoncross,:breedcross]) && push!(popidls,popid)
  end
  issubset(popidls,ped.member) || @error string("inconsistent popid=",popidls)
  filter!(row->in(row[:member],popidls), magicped.offspringinfo)
  magicped.designinfo = ped 
  magicped
end

function magicped2ped(magicped::MagicPed; isfounderinbred::Bool=true)
    if isa(magicped.designinfo, Dict{String, DesignInfo})                
        infols = [[key,join(val.founders),val.designtype,val.pedigree] for (key,val) in magicped.designinfo]
        b = [i[3] == :juncdist for i in infols]
        infols = infols[.!b]
        isempty(infols) && return nothing
        o = sortperm([i[2] for i in infols])
        infols = infols[o]
        b = isnothing.(last.(infols))
        if any(b)
            msg = string("unknown pedigrees for subpopulations: ", (first.(infols))[b])
            @warn msg
        end        
        pedigrees = Vector{Pedigree{String}}((last.(infols))[.!b])
        ped = Pedigrees.mergeped(pedigrees; isfounderinbred)
    elseif isa(magicped.designinfo, Pedigree)
        ped = magicped.designinfo
    else
        @error string("Could not extract ped from designinfo=",magicped.designinfo)
        return nothing
    end
    ped
end

function split_subpop(magicpedfile::AbstractString;
    outstem::AbstractString = "outstem",
    workdir::AbstractString=pwd())
    magicped = readmagicped(magicpedfile)
    split_subpop(magicped; outstem, workdir)
end

function split_subpop(magicped::MagicPed;
    outstem::AbstractString = "outstem",
    workdir::AbstractString=pwd())
    memls = unique(magicped.offspringinfo[!,:member])
    res = String[]
    for mem in memls
        subped = submagicped(magicped,mem)
        mem2 = replace(mem,"/"=>"_","\\"=>"_")
        outfile = getabsfile(workdir,outstem*"_"*mem2*"_ped.csv")
        savemagicped(outfile, subped)
        push!(res,outfile)
    end
    res
end

function rename_nonfounder!(magicped::MagicPed; 
    suffix::AbstractString=Random.randstring(),
    keeppopid::Bool=true)
    keepmembers = keeppopid ? unique(magicped.offspringinfo[!,:member]) : []
    magicped.designinfo = Pedigrees.rename_nonfounder(magicped.designinfo;suffix,keepmembers)
    if !keeppopid    
        magicped.offspringinfo[!,:member] .= [string(i,"_",suffix) for i in magicped.offspringinfo[!,:member]]
    end
    magicped
end

function pad_offspringdf!(offspringdf::DataFrame)
    cols = names(offspringdf)
    issubset(["individual","member"],cols) || @error string("member and designcode columns must be in offspringdf; cols=",cols)
    in("ishomozygous",cols) || insertcols!(offspringdf,3,:ishomozygous => false)
    in("isfglexch",cols) || insertcols!(offspringdf,4,:isfglexch => false)
    offspringdf
end

function pad_designdf!(designdf::DataFrame)
    cols = names(designdf)
    issubset(["member","designcode"],cols) || @error string("member and designcode columns must be in designdf; cols=",cols)
    in("founders",cols) || insertcols!(designdf,3,:founders => "NA")    
    designdf
end

"""
    merge_pedfiles(pedfiles; keyargs...)

merge pedfiles into a single pedfile. 

# Positional arguments

`pedfiles::AbstractString`: a list of pedfiles. 

# Keyword arguments

`isped::Bool`: if true, designinfo in pedfiles is in form of pedigree, and otherwise designcode. 

`isfounderinbred::Bool=true`: if true, founders are inbred. 

`outstem::AbstractString=popid`: stem of output filename.

`workdir::AbstractString=pwd()`: directory for reading and saving files.

"""
function merge_pedfiles(pedfiles::AbstractVector;     
    isped::Bool,
    isfounderinbred::Bool=true,
    workdir::AbstractString=pwd(),
    outstem::AbstractString="outstem")     
    if isped
        merge_pedfiles_ped(pedfiles; isfounderinbred, workdir,outstem)        
    else
        merge_pedfiles_designcode(pedfiles; workdir,outstem)        
    end

end
function merge_pedfiles_designcode(pedfiles::AbstractVector;     
    commentstring::AbstractString="##",
    workdir::AbstractString=pwd(),
    outstem::AbstractString="outstem")        
    pedfilels = [getabsfile(workdir,i) for i in pedfiles]
    b = isfile.(pedfilels)
    all(b) || @error string(pedfiles[.!b], " do not exist!")
    delim=','
    peddict = MagicBase.readmultitable(pedfiles[1]; delim,commentstring,isparse=false)
    if !issubset(["designinfo","offspringinfo"],keys(peddict))
        error(string("keys designinfo and offspringinfo are not in ",keys(peddict)))
    end    
    designdf = pad_designdf!(peddict["designinfo"])
    offdf = pad_offspringdf!(peddict["offspringinfo"])
    popfirst!(pedfilels)
    for pedfile in pedfilels
        @info string("mering pedfile=",pedfile)
        peddict = MagicBase.readmultitable(pedfile; delim,commentstring,isparse=false)
        if !issubset(["designinfo","offspringinfo"],keys(peddict))
            error(string("keys designinfo and offspringinfo are not in ",keys(peddict)))
        end    
        designdf2 = pad_designdf!(peddict["designinfo"])
        offdf2 = pad_offspringdf!(peddict["offspringinfo"])
        # designdf = vcat(designdf,designdf2)
        # offdf = vcat(offdf,offdf2)
        append!(designdf, designdf2)
        append!(offdf, offdf2)
    end
    dup = filter(x->x.second>1, countmap(designdf[!,:member]))
    isempty(dup) || @error string("duplicated design: ",dup)
    dup = filter(x->x.second>1, countmap(offdf[!,:individual]))
    isempty(dup) || @error string("duplicated offspring: ",dup)
    if unique(designdf[!,:founders]) == ["NA"]
        select!(designdf,Not([:founders]))
    end
    if unique(offdf[!,:isfglexch]) == [false]
        select!(offdf,Not([:isfglexch]))
        if unique(offdf[!,:ishomozygous]) == [false]
            select!(offdf,Not([:ishomozygous]))
        end
    end
    pedfile = getabsfile(workdir,outstem*"_ped.csv")
    open(pedfile,"w") do io
        write(io,"RABBIT,designinfo\n")        
        CSV.write(io, designdf; append=true,writeheader=true)
        write(io,"RABBIT,offspringinfo\n")    
        CSV.write(io, offdf; append=true,writeheader=true)
    end
    pedfile
end

function merge_pedfiles_ped(pedfiles::AbstractVector; 
    isfounderinbred::Bool=true,
    workdir::AbstractString=pwd(),
    outstem::AbstractString="outstem")
    pedfilels = [getabsfile(workdir,i) for i in pedfiles]
    b = isfile.(pedfilels)
    all(b) || @error string(pedfiles[.!b], " do not exist!")
    magicped = readmagicped(pedfilels[1])
    rename_nonfounder!(magicped; keeppopid=true,suffix="1")
    peddf = Pedigrees.ped2df(magicped.designinfo)
    finfo = magicped.founderinfo
    offinfo = magicped.offspringinfo
    popfirst!(pedfilels)
    for (i,pedfile) in enumerate(pedfilels)
        magicped2 = readmagicped(pedfile)
        rename_nonfounder!(magicped2;keeppopid=true,suffix=string(i+1))
        peddf2 = Pedigrees.ped2df(magicped2.designinfo)
        peddf = unique(vcat(peddf,peddf2))
        finfo = unique(vcat(finfo,magicped2.founderinfo))
        offinfo = unique(vcat(offinfo,magicped2.offspringinfo))
    end
    sort!(peddf,:generation)
    if isfounderinbred
        f1df = peddf[peddf[!,:generation] .== 1,:]
        dict = Dict()
        for row in eachrow(f1df)
            push!(dict,row[:member]=>join(sort(Vector(row[[:mother,:father]])),"/"))
        end
        peddf[!,[:member,:mother,:father]] .= [get(dict,i,i) for i in Matrix(peddf[!,[:member,:mother,:father]])]
    end
    unique!(peddf)    
    ped = Pedigree(peddf)
    if ped.member[1:ped.nfounder] != finfo[!,:individual]
        @error string("inconsistent founderid between ped and founderinfo")
    end
    memls = unique(offinfo[!,:member])
    d = setdiff!(memls,ped.member)
    isempty(d) || @error string("popid=",d, " not in ped")
    magicped = MagicPed(ped,finfo,offinfo)
    savemagicped(getabsfile(workdir,outstem*"_ped.csv"),magicped)
end


function get_founder2offspring(magicped::MagicPed; isindex::Bool=false)
    subpop_offspring = get_subpop2offspring(magicped; isindex)
    founder_subpop = get_founder2subpop(magicped)
    Dict([p=>reduce(vcat,[subpop_offspring[i] for i in subpopls])
        for (p,subpopls) in founder_subpop])
end

function get_subpop2offspring(magicped::MagicPed;
    isindex::Bool=false)
    if isindex
        ls = magicped.offspringinfo[!,:member]
        Dict([i => findall(ls .== i) for i in unique(ls)])
    else
        offls = magicped.offspringinfo[!,:individual]
        ls = magicped.offspringinfo[!,:member]
        Dict([i => offls[ls .== i] for i in unique(ls)])
    end
end

function get_founder2subpop(magicped::MagicPed)
    subpop_founder = get_subpop2founder(magicped)
    dict = Dict{String,Vector{String}}()
    for (subpop,pp) in subpop_founder
        for p in pp
            if haskey(dict,p)
                push!(dict[p] ,subpop)
            else
                push!(dict, p=>[subpop])
            end
        end
    end
    dict
end

function get_subpop2founder(magicped::MagicPed;isindex::Bool=false)
    if isindex 
        founder2index = Dict(magicped.founderinfo[!,:individual] .=> 1:size(magicped.founderinfo,1))
    end
    design = magicped.designinfo
    memls = unique(magicped.offspringinfo[!,:member])
    if isa(design,Pedigree)         
        Dict([begin
            subped = getsubped(design,convert(String,i))
            fls = subped.member[1:subped.nfounder]
            i => (isindex ? [founder2index[f] for f in fls] : fls)
        end  for i in memls])
    elseif isa(design,Dict{String,DesignInfo})         
        Dict([convert(String,popid)  => (isindex ? [founder2index[f] for f in subdesign.founders]  : subdesign.founders) for (popid,subdesign) in design])
    else
        @error string("unknown designinfo type: ",typeof(design))
    end    
end

function generate_magicped(;
    designcodes::AbstractVector = ["P1/P2=>DH","P4/3/P2/P3//P5/P6=>3"],
    founders::AbstractVector = [missing for _ in 1:length(designcodes)],
    subpopsizes::AbstractVector =10*ones(Int, length(designcodes)),
    workdir::AbstractString=pwd(),
    outstem::Union{Nothing,AbstractString}=nothing)
    npop = length(designcodes) 
    npop == length(founders) || @error string("inconsistent #founders between designcodes and founders")
    npop == length(subpopsizes) || @error string("inconsistent #founders between designcodes and subpopsizes")
    members = string.("pop",1:npop)
    founders2 = [i in ["missing","NA"] ? "NA" : i for i in string.(founders)]
    designdf = DataFrame(member=members,founders=founders2,designcode=designcodes)        
    offdf = DataFrame(individual=String[],member=String[])
    for i in eachindex(members)
        subsize  = round(Int, subpopsizes[i])
        df = DataFrame(individual = string.(members[i],"_", 1:subsize), member = members[i])
        offdf = vcat(offdf,df)
    end    
    peddict =Dict("designinfo"=>designdf,"offspringinfo"=>offdf)
    magicped = MagicBase.parsemagicped(peddict)
    if !isnothing(outstem)
        try 
            plotmagicped(magicped)
            savefig(getabsfile(workdir, outstem*"_ped.png"))
        catch err
            msg = string(err,". Could not plot magicped")
            @warn msg            
        end
        outfile = getabsfile(workdir, outstem*"_ped.csv")
       savemagicped(outfile, magicped)
    end
    magicped
end


function get_connected_pops(magicped::MagicPed)
    founderls = magicped.founderinfo[!,:individual]
    cc = get_connected_parents(magicped)
    ccparentls = [founderls[i] for i in cc]
    f2pop = MagicBase.get_founder2subpop(magicped)
    ccpopls = [unique(reduce(vcat,[f2pop[i] for i in c])) for c in ccparentls]
    popls = unique(magicped.offspringinfo[!,:member])
    sum(length.(ccpopls)) == length(popls) || error("unexpected connected subpopulations!")
    issetequal(popls, reduce(vcat,ccpopls)) || error("unexpected connected subpopulations!")
    ccpopls
end


function get_connected_parents(magicped::MagicPed)
    pop2f = MagicBase.get_subpop2founder(magicped; isindex = true)
    np = size(magicped.founderinfo,1)
    adj = zeros(Int, np,np)
    for (subpop,pp) in pop2f
        adj[pp,pp] .= 1
    end
    adj=sign.(adj + adj')
    g = SimpleGraph(adj)
    connected_components(g)
end

