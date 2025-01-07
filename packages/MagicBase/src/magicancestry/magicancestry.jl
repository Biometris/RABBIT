"""
    MagicAncestry

mutable struct that stores the results of haplotype reconstruction by magicreconstruct. 

    MagicAncestry(magicped,markermap,foundergeno,statespace,
        viterbipath, diploprob,genoprob,haploprob,loglike,misc)

inner constructor. See also [`readmagicancestry`](@ref).

# Fields

`magicped::MagicPed`: breeding pedigree information. See also [`MagicPed`](@ref).

`markermap::Vector{DataFrame}`: marker map for each chromosome. `markermap[c]`
gives the marker map of chromosome c. 

`foundergeno::Vector{Matrix}`: haplotypes in founders. `foundergeno[c][m,f]`
gives the founder f at marker m of chromosome c.

`statespace::AbstractDict`: definition of ancestral haplotype states corresponding to `haploprob`.

`viterbipath::Union{Nothing,Vector{Matrix}}`: optimal state path obtained by the Viterbi
algorithm. `viterbipath[c][m,o]` gives the ancestral diplotype state of offspring o at marker m for chromosome c.

`diploprob::Union{Nothing,Vector{Vector{SparseArrays.SparseMatrixCSC}}}`:
marginal posterior probabilities for diplotypes. `diploprob[c][o][m,s]` gives the
probability of offspring o at ancestral diplotype state s for marker m of chromosome c.

`genoprob::Union{Nothing,Vector{Vector{SparseArrays.SparseMatrixCSC}}}`:
marginal posterior probabilities for genotypes. `genoprob[c][o][m,s]` gives the
probability of offspring o at ancestral genotype state s for marker m of chromosome c.

`haploprob::Union{Nothing,Vector{Vector{SparseArrays.SparseMatrixCSC}}}`:
marginal posterior probabilities for haplotypes. `genoprob[c][o][m,s]` gives the
probability of offspring o at ancestral haplotype state s for marker m of chromosome c.

`inbredcoef::Union{Nothing,Vector{Matrix}}`: realized inbreeding coefficients. 
`inbredcoef[c][m,o]` gives the inbreeding coefficients at marker m of chromosome c in offpsring o. 

`loglike::AbstractMatrix`: loglike[c,o] gives the log likelihood for chromsome c
of offspring o.

`misc::Dict{String, DataFrame}`: contains information such as (1) "deletion" dataframe
for markers that were removed from markermap, (2) "correction" dataframe for
parental error correction.

"""
mutable struct MagicAncestry
    magicped::MagicPed
    # markermap[:markerid,:chrid,:poscm]
    markermap::Vector{DataFrame}
    # foundergeno[chr][marker,founder]
    foundergeno::AbstractVector
    statespace::AbstractDict    
    # viterbipath[chr][:,offspring]: most probably path for offspring at chr.
    viterbipath::Union{Nothing,Vector{Matrix}}
    # diploprob[chr][offspring][marker, diplostate]
    diploprob::Union{Nothing,Vector{Vector{SparseArrays.SparseMatrixCSC}}}
    # genoprob[chr][offspring][marker, genostate]
    genoprob::Union{Nothing,Vector{Vector{SparseArrays.SparseMatrixCSC}}}
    # haploprob[chr][offspring][marker, haplostate]
    haploprob::Union{Nothing,Vector{Vector{SparseArrays.SparseMatrixCSC}}}
    # inbredcoef[chr][marker, offspring]
    inbredcoef::Union{Nothing,Vector{Matrix}}
    # loglike[chr,offspring]
    loglike::AbstractMatrix
    misc::Dict{String, DataFrame}
    function MagicAncestry(magicped,markermap,foundergeno,statespace,
        viterbipath,diploprob,genoprob,haploprob,inbredcoef, loglike,misc)
        dims = [length(markermap),size(foundergeno,1),size(loglike,1)]
        if length(union(dims))!=1
            error("inconsistent number of linkage groups ",dims)
        end
        lgsize = [size(i,1) for i=markermap]
        dims = [lgsize,size.(foundergeno,1)]
        if length(union(dims))!=1
            error("inconsistent number of markers within linkage groups ",dims)
        end
        dims = [size(magicped.offspringinfo,1),size(loglike,2)]
        if length(union(dims))!=1
            error("inconsistent number of offspring ",dims)
        end
        nfhaplo=size(first(foundergeno),2)
        nfounder=size(magicped.founderinfo,1)
        if !(nfounder==nfhaplo || 2*nfounder==nfhaplo)
            error("inconsistent number of founders")
        end
        new(magicped,markermap,foundergeno,statespace,
            viterbipath,diploprob,genoprob,haploprob,
            inbredcoef, loglike,misc)
    end
end

"""
    savemagicancestry(outputfile,magicancestry,workdir=pwd())

save magicancestry into outputfile.

# Positional arguments

`outputfile::AbstractString`: file for saving magicancestry. 

`magicancestry::MagicAncestry`: a struct returned by magicreconstruct. 

# Keyward arguments

`workdir::AbstractString = pwd()`: directory for writing outputfile.

"""
function savemagicancestry(outputfile::AbstractString,magicancestry::MagicAncestry;
    workdir::AbstractString = pwd(),
    tempdirectory::AbstractString = tempdir())
    ext = last(split_allext(outputfile))
    in(ext, [".jld2.tar.gz", ".csv.gz.tar", ".csv.gz",".csv"]) || @error string("ext=",
        ext, ", not .csv.gz, .csv, or .csv.gz.tar")
    if in(ext, [".csv",".csv.gz"])
        savemagicancestry_csvgz(outputfile,magicancestry; workdir)
    elseif ext == ".csv.gz.tar"
        savemagicancestry_csvgztar(outputfile,magicancestry;
            workdir,tempdirectory)
    elseif ext == ".jld2.tar.gz"
        savemagicancestry_jld2targz(outputfile,magicancestry;
            workdir,tempdirectory)
    else
        error(string("unkown fileext=",ext))
    end
end

function get_ancestry_chrlab(chr::Integer, nchr::Integer)
    strlen = ndigits(nchr)
    "chr"*lpad(string(chr),strlen,"0")
end

function savemagicancestry_jld2targz(outputfile::AbstractString,magicancestry::MagicAncestry;
    workdir::AbstractString = pwd(),
    tempdirectory::AbstractString = tempdir())
    ext = last(split_allext(outputfile))
    in(ext, [".jld2.tar.gz", ".tar.gz", ".tar"]) || @error string("ext=",
        ext, ", not in [.jld2.tar.gz, .tar.gz, .tar]")
    nchr = length(magicancestry.markermap)
    jltempdir = mktempdir(tempdirectory;
        prefix="jl_savemagicancestry_", cleanup=true)
    try        
        outstem = first(split_allext(basename(outputfile)))
        for chr in 1:nchr
            chrancestry = MagicBase.submagicancestry(magicancestry;chrsubset=[chr])
            chrlab = get_ancestry_chrlab(chr,nchr)
            chroutfile = joinpath(jltempdir,string(outstem,"_",chrlab, ".jld2"))            
            chrid = magicancestry.markermap[chr][1,:linkagegroup]		
            JLD2.save(chroutfile, chrid, chrancestry)
        end
        outputfile2 = getabsfile(workdir, outputfile)
        if ext == ".tar"
            create_tar(jltempdir, outputfile2)
        else
            create_targz(jltempdir, outputfile2)
        end
        outputfile2
    finally
        rm.(readdir(jltempdir;join=true);force=true)
        rm(jltempdir; force=true,recursive=true)
    end
end

function savemagicancestry_csvgztar(outputfile::AbstractString,magicancestry::MagicAncestry;
    workdir::AbstractString = pwd(),
    tempdirectory::AbstractString = tempdir())
    ext = last(split_allext(outputfile))
    in(ext, [".tar"]) || @error string("ext=",ext, ", not .tar]")
    nchr = length(magicancestry.markermap)
    jltempdir = mktempdir(tempdirectory;
        prefix="jl_savemagicancestry_", cleanup=true)
    try        
        outstem = first(split_allext(basename(outputfile)))
        for chr in 1:nchr
            chrancestry = MagicBase.submagicancestry(magicancestry;chrsubset=[chr])
            chrlab = get_ancestry_chrlab(chr,nchr)
            chroutfile = joinpath(jltempdir,string(outstem,"_",chrlab, ".csv.gz"))
            savemagicancestry(chroutfile,chrancestry)            
        end
        outputfile2 = getabsfile(workdir, outputfile)
        create_tar(jltempdir, outputfile2)
        outputfile2
    finally
        rm.(readdir(jltempdir;join=true);force=true)
        rm(jltempdir; force=true,recursive=true)
    end
end

function savemagicancestry_csvgz(outputfile::AbstractString,magicancestry::MagicAncestry;
    workdir::AbstractString = pwd())
    delim = ','
    ext = last(split_allext(outputfile))
    in(ext, [".csv.gz",".csv"]) || @error string("ext=",
        ext, ", not in [.csv.gz, .csv]")
    myopen = ext == ".csv.gz" ? GZip.open : open
    outputfile2 = getabsfile(workdir,outputfile)
    myopen(outputfile2, "w") do io
        initial = "RABBIT"
        ped=magicancestry.magicped.designinfo
        df = MagicBase.designinfo2df(ped)        
        appenddf(io, df; delim,initial,dfname="designinfo")
        appenddf(io, magicancestry.magicped.founderinfo; delim,initial,dfname="founderinfo")
        appenddf(io, magicancestry.magicped.offspringinfo; delim,initial,dfname="offspringinfo")
        # appenddf(io, magicancestry.deletion; delim,initial,dfname="deletion")
        # appenddf(io, magicancestry.correction; delim,initial,dfname="correction")
        savefoundergeno(io, magicancestry; delim, initial)
        savestatespace(io, magicancestry; delim, initial)
        savecondprob_ind(io, magicancestry; delim, initial,probtype="diploprob")
        savecondprob_ind(io, magicancestry; delim, initial,probtype="genoprob")
        savecondprob_ind(io, magicancestry; delim, initial,probtype="haploprob")
        saveinbredcoef(io, magicancestry; delim, initial)
        saveviterbipath(io, magicancestry; delim, initial)
        saveloglike(io, magicancestry; delim, initial)
    end
    outputfile2
end

function appenddf(io::IO, df::AbstractDataFrame;
    delim::AbstractChar=',',
    dfname::AbstractString="dfname",
    initial::AbstractString="RABBIT",
    header::Bool=true)
    write(io,initial,delim,dfname,"\n")
    CSV.write(io, df; delim,header,append=true)
end

function savefoundergeno(io::IO, magicancestry::MagicAncestry;
    delim::AbstractChar,
    initial::AbstractString="initial")
    nchr = length(magicancestry.foundergeno)
    fkind = only(unique(reduce(vcat,[unique(i[!,:founderformat]) for i=magicancestry.markermap])))
    if fkind == "GT_haplo"
        fhaplo=[magicancestry.foundergeno[ch] for ch=1:nchr]
    elseif fkind == "GT_phased"
        fhaplo=[join.(magicancestry.foundergeno[ch],"|") for ch=1:nchr]
    else
        @error string("unknown founder genotype format: ",fkind)
    end
    write(io,initial,delim,"foundergeno","\n")
    rowid = [names(magicancestry.markermap[1]); magicancestry.magicped.founderinfo[!,:individual]]
    join(io, rowid, delim)
    write(io,"\n")
    for chr=1:nchr
        chrmap = Matrix(magicancestry.markermap[chr])
        fhaplo2=hcat(chrmap,fhaplo[chr])
        writedlm(io,fhaplo2,delim)
    end
    flush(io)
end

function savestatespace(io::IO, magicancestry::MagicAncestry;
    delim::AbstractChar,
    initial::AbstractString="initial")
    states = magicancestry.statespace
    fgls = states["haplotype"]
    nfgl = length(fgls)
    haplodf = DataFrame(hapotypeindex=1:length(fgls), haploindex=states["haploindex"], haplotype=fgls)   
    appenddf(io, haplodf; delim,initial,dfname=names(haplodf)[3])    
    if isnothing(magicancestry.viterbipath)
        if isnothing(magicancestry.diploprob) && isnothing(magicancestry.haploprob) 
            @error string("Could not detect if the state space is the set of diplotypes")
        else
            isdiplo = !isnothing(magicancestry.diploprob) 
        end
    else
        isdiplo = maximum([maximum(i) for i in magicancestry.viterbipath]) == nfgl^2
    end    
    if isdiplo 
        genoindex = [[i,j] for i=1:nfgl for j=1:i]
        genotype=[join(fgls[i],"|") for i in genoindex]
        genodf = DataFrame(genotypeindex=1:length(genotype), genoindex=[join(i,"|") for i in genoindex], genotype=genotype)
        diploindex = [[i,j] for i=1:nfgl for j=1:nfgl]
        diplotype=[join(fgls[i],"|") for i in diploindex]
        diplodf = DataFrame(diplotypeindex=1:length(diplotype), diploindex=[join(i,"|") for i in diploindex], diplotype=diplotype)
        for df=[genodf,diplodf]
            appenddf(io, df; delim,initial,dfname=names(df)[3])
        end
    else
        write(io,initial,delim,"genotype","\n")
        write(io, "nothing\n")
        write(io,initial,delim,"diplotype","\n")
        write(io, "nothing\n")
    end
end

function savecondprob_ind(io::IO, magicancestry::MagicAncestry;
    probtype::AbstractString="haploprob",
    delim::AbstractChar=',',
    initial::Union{Nothing,AbstractString}=nothing)
    if probtype=="haploprob"
        condprob = magicancestry.haploprob
        Jlabel = "haplotypeindex"
        nJlabel = "nhaplotype"
    elseif probtype == "genoprob"
        condprob = magicancestry.genoprob
        Jlabel = "genotypeindex"
        nJlabel = "ngenotype"
    elseif probtype=="diploprob"
        condprob = magicancestry.diploprob
        Jlabel = "diplotypeindex"
        nJlabel = "ndiplotype"
    else
        error(string("unknown probtype:",probtype))
    end
    isnothing(initial) || write(io,initial,delim,probtype,"\n")
    if isnothing(condprob)
        write(io, "nothing\n")
        return 0
    end
    offls = magicancestry.magicped.offspringinfo[!,:individual]
    noff = length(offls)    
    nchr = length(magicancestry.markermap)
    rowid =  ["linkagegroup","offspring", "nmarker", nJlabel, "markerindex",Jlabel,probtype]
    join(io, rowid, delim)
    write(io,"\n")    
    nthread = Threads.nthreads()
    offblockls = Iterators.partition(1:noff, 3*nthread)
    for chr=1:nchr
        chrid = magicancestry.markermap[chr][1,:linkagegroup]
        chrcondprob = condprob[chr]
        nsnp, nstate = size(chrcondprob[1])        
        for offblock in offblockls
            res = Vector{String}(undef,length(offblock))
            ThreadsX.foreach(eachindex(offblock)) do i
                off = offblock[i]
                if issparse(chrcondprob[off])
                    I, J, V=findnz(chrcondprob[off])
                else
                    sparseprob = sparse(round.(chrcondprob[off],digits=4))
                    I, J, V=findnz(sparseprob)                     
                end
                length(I) == length(J) == length(V) || @error string("inconsistent lengths of I, J, V")
                ls = [chrid, offls[off], nsnp, nstate,join(I, "|"), join(J, "|"), join(V, "|")]                            
                res[i] = join(ls,delim)
            end
            for i in res
                write(io,i, "\n")
                flush(io)
            end
        end
    end
    flush(io)
end


function tostringprob(prob::SparseVector)
    I, V =  findnz(prob)
    join(I,"|")*"=>"*join(V,"|")
end

function tostringprob(prob::SparseMatrixCSC)
    [tostringprob(prob[i,:]) for i=1:size(prob,1)]
end

function savecondprob_marker(io::IO, magicancestry::MagicAncestry;
    probtype::AbstractString="haploprob",
    delim::AbstractChar=',',
    initial::Union{Nothing,AbstractString}=nothing)
    if probtype=="haploprob"
        condprob = magicancestry.haploprob
    elseif probtype == "genoprob"
        condprob = magicancestry.genoprob
    elseif probtype=="diploprob"
        condprob = magicancestry.diploprob
    else
        error(string("unknown probtype:",probtype))
    end
    isnothing(initial) || write(io,initial,delim,probtype,"\n")
    if isnothing(condprob)
        write(io, "nothing\n")
        return 0
    end
    offls = magicancestry.offspringinfo[!,:individual]
    rowid = [names(magicancestry.markermap[1]); offls]
    join(io, rowid, delim)
    write(io,"\n")
    for chr=1:length(magicancestry.markermap)
        chrmap = Matrix(magicancestry.markermap[chr])        
        res = reduce(hcat,tostringprob.(condprob[chr]))
        res2=hcat(chrmap,res)
        writedlm(io,res2,delim)
    end
    flush(io)
end

function saveinbredcoef(io::IO,magicancestry::MagicAncestry;
    delim::AbstractChar=',',
    initial::AbstractString="initial")    
    if isnothing(magicancestry.inbredcoef)
        write(io,initial,delim,"inbredcoef","\n")
        write(io, "nothing\n")
        return 0
    end
    write(io,initial,delim,"inbredcoef","\n")
    indls = Symbol.(magicancestry.magicped.offspringinfo[!,:individual])
    nchr = length(magicancestry.markermap)
    for chr in 1:nchr
        chrid = magicancestry.markermap[chr][1,:linkagegroup]
        if !isnothing(magicancestry.viterbipath)
            inbredcoef = [join(Int.(i),"|") for i in eachcol(magicancestry.inbredcoef[chr])]
        else
            inbredcoef = [join(i,"|") for i in eachcol(magicancestry.inbredcoef[chr])]
        end
        df = DataFrame(chromosome=chrid,offspring=indls,inbredcoef=inbredcoef)
        CSV.write(io, df; delim,header = chr==1,append=true)
    end
    nothing
end

function saveviterbipath(io::IO,magicancestry::MagicAncestry;
    delim::AbstractChar=',',
    initial::AbstractString="initial")    
    if isnothing(magicancestry.viterbipath)
        write(io,initial,delim,"viterbipath","\n")
        write(io, "nothing\n")
        return 0
    end
    # res = [[HMM.tostringpath(HMM.tojumppath(i)) for i=eachcol(j)] for j=magicancestry.viterbipath]
    res = [[tostringpath(tojumppath(i)) for i=eachcol(j)] for j=magicancestry.viterbipath]    
    res2=[join(i,"|") for i=eachrow(reduce(hcat,res))]
    offls = magicancestry.magicped.offspringinfo[!,:individual]
    df = DataFrame([offls res2],[:individual,:viterbipath])
    appenddf(io, df; delim,initial,dfname="viterbipath")
end

function saveloglike(io::IO,magicancestry::MagicAncestry;
    delim::AbstractChar=',',
    initial::AbstractString="initial")
    indls = magicancestry.magicped.offspringinfo[!,:individual]
    df= DataFrame(magicancestry.loglike,Symbol.(indls))
    chridls = [i[1,:linkagegroup] for i= magicancestry.markermap]
    df2 = hcat(DataFrame(:linkagegroup=>chridls),df)
    appenddf(io, df2; delim,initial,dfname="loglike")
end

"""

    readmagicancestry(ancestryfile, workdir=pwd(),tempdirectory=tempdir())

read ancestryfile in the directory workdir and return magicancestry::MagicAncestry. See [`MagicAncestry`](@ref). 

# Positional argument

`ancestryfile`: file storing magicancestry that is generated by [`savemagicancestry`](@ref). It results from `magicreconstruct`.

# Keyward arguments

`workdir::AbstractString = pwd()`: directory for reading ancestryfile,

`tempdirectory::AbstractString=tempdir()`: temparary directory. 

"""
function readmagicancestry(ancestryfile::AbstractString;
    workdir::AbstractString=pwd(),
    tempdirectory::AbstractString=tempdir())
    ext = last(split_allext(ancestryfile))    
    if in(ext, [".csv",".csv.gz"])
        readmagicancestry_csvgz(ancestryfile;workdir)
    elseif ext == ".csv.gz.tar"
        readmagicancestry_csvgztar(ancestryfile;workdir)
    elseif ext == ".jld2.tar.gz"
        readmagicancestry_jld2targz(ancestryfile;workdir,tempdirectory)
    else
        error(string("unknow fileext=",ext))
    end
end

function readmagicancestry_jld2targz(ancestryfile::AbstractString;
    workdir::AbstractString=pwd(),
    tempdirectory::AbstractString=tempdir())
    jltempdir = mktempdir(tempdirectory;
        prefix="jl_readmagicancestry_", cleanup=true)
    try
        extract_targz(getabsfile(workdir,ancestryfile),jltempdir)
        filels = readdir(jltempdir;sort=true,join=true)
        ancestryls = [only(JLD2.load(file)).second for file in filels]
        MagicBase.mergeancestry(ancestryls)
    finally
        rm.(readdir(jltempdir;join=true);force=true)
        rm(jltempdir; force=true,recursive=true)
    end
end

function readmagicancestry_csvgztar(ancestryfile::AbstractString;
    workdir::AbstractString=pwd(),
    tempdirectory::AbstractString=tempdir())
    jltempdir = mktempdir(tempdirectory;
        prefix="jl_readmagicancestry_", cleanup=true)
    try
        extract_tar(getabsfile(workdir,ancestryfile),jltempdir)
        filels = readdir(jltempdir;sort=true,join=true)
        ancestryls = [readmagicancestry(file) for file in filels]
        MagicBase.mergeancestry(ancestryls)
    finally
        rm.(readdir(jltempdir;join=true);force=true)
        rm(jltempdir; force=true,recursive=true)
    end
end

function readmagicancestry_csvgz(ancestryfile::AbstractString;    
    commentstring::AbstractString="##",
    workdir::AbstractString=pwd())
    delim=','    
    ancestryfile2 = getabsfile(workdir,ancestryfile)
    isfile(ancestryfile2) || @error string(ancestryfile2, " does not exist")
    # issparse = true, CSV.read is used to parse string, when cells with too long string will be truncated
    res= readmultitable(ancestryfile2;delim,commentstring,isparse=false)
    statespace= parsestatespace(res)
    magicped = parsemagicped(res)
    misc = Dict{String, DataFrame}()
    for (key,val) in res
        key[1:5]=="misc|" && push!(misc, key[6:end]=>val)
    end
    genodf = res["foundergeno"]
    parsecol_pos_error!(genodf)        
    founderid = magicped.founderinfo[!,:individual]
    isfounderinbred = length(founderid) == length(statespace["haplotype"])    
    # newfounderid  are assumed to be consistent with founderid and offspringid
    mapdf, fgenodf, _, offgenodf, _ = splitgenodf(genodf,founderid,[];
        isphysmap=false)    
    markermap, foundergeno, _ = parsegenodf!(mapdf, fgenodf, offgenodf; isfounderinbred)    
    empty!(genodf)
    viterbidf = res["viterbipath"]
    viterbipath = isempty(viterbidf) ? nothing : parseviterbi(viterbidf)
    diploprobdf = res["diploprob"]
    diploprob = isempty(diploprobdf) ? nothing : parsecondprob_ind!(diploprobdf)
    genoprobdf = res["genoprob"]
    # nfgl*(nfgl+1) ÷ 2
    genoprob = isempty(genoprobdf) ? nothing : parsecondprob_ind!(genoprobdf)
    haploprobdf = res["haploprob"]
    haploprob = isempty(haploprobdf) ? nothing : parsecondprob_ind!(haploprobdf)    
    if haskey(res,"inbredcoef")
        inbredcoef = isempty(res["inbredcoef"]) ? nothing : parse_inbredcoef(res["inbredcoef"])
    else
        inbredcoef = nothing
    end
    loglike=Matrix(res["loglike"][!,2:end])
    MagicAncestry(magicped,markermap,foundergeno,statespace,
        viterbipath,diploprob,genoprob,haploprob,inbredcoef,loglike,misc)
end

function parseviterbi(viterbidf::DataFrame)
    a=split.(viterbidf[!,:viterbipath],"|")
    res = [[todiscretepath(tojumppath(i)) for i=j] for j=a]
    res2 = reduce(hcat, res)
    [reduce(hcat, i) for i=eachrow(res2)]    
end

function parsecondprob_ind!(condprobdf::DataFrame)
    # rowid =  ["linkagegroup","offspring","nmarker", "nstate",
    #           "markerindex","stateindex","probtype"]
    for j in 3:4
        if eltype(condprobdf[!,j]) <: AbstractString
            condprobdf[!,j] .= parse.(Int,condprobdf[!,j])
        end
    end
    groupdf = groupby(condprobdf,1)
    [begin
        mtx=Matrix(groupdf[ch][!,3:7])
        for col in 3:4
            mtx[:,col] .= [parse.(Int,split(i,"|")) for i = mtx[:,col]]
        end
        col=5
        mtx[:,col] .= [parse.(Float64,split(i,"|")) for i = mtx[:,col]]
        [sparse(i[3],i[4],i[5],i[1],i[2]) for i=eachrow(mtx)]
    end for ch=1:length(groupdf)]
end

function parse_inbredcoef(inbredcoef::DataFrame)
    # rowid =  ["linkagegroup","offspring","inbredcoef"]
    isbool = issubset(unique(split(inbredcoef[1,:inbredcoef],"|")),["0","1"])
    parsetype = isbool ? Bool : Float64
    groupdf = groupby(inbredcoef,1)
    [reduce(hcat,[parse.(parsetype,split(i,"|")) for i in df[!,:inbredcoef]]) for df in groupdf]
end

function parseprobcell(x::AbstractString,nstate::Integer)
    ls=split.(split(x,"=>"),"|")
    sparsevec(parse.(Int,ls[1]),parse.(Float64,ls[2]),nstate)
end

function parsecondprob_marker(condprobdf::DataFrame,nstate::Integer)
    groupdf = groupby(condprobdf,2)
    [begin
        mtx=Matrix(groupdf[ch][:,4:end])
        mtx2=parseprobcell.(mtx,nstate)
        [sparse(reduce(hcat,col)') for col=eachcol(mtx2)]
    end for ch=1:length(groupdf)]
end

function parsestatespace(res::AbstractDict)
    haploindex = res["haplotype"][!,:haploindex]
    haploindex2 = eltype(haploindex) <: AbstractString ? parse.(Int,haploindex) : haploindex
    haplotype = string.(strip.(res["haplotype"][!,:haplotype]))    
    statespace = Dict("haplotype"=>haplotype,"haploindex"=>haploindex2)
    statespace
end


"""
    setgenoprob!(magicancestry::MagicAncestry)

set magicancestry.genoprob from magicancestry.diploprob.
"""
function setgenoprob!(magicancestry::MagicAncestry; posteriordigits)
    if isnothing(magicancestry.diploprob)
        error("diploprob is required")
    end
    statespace=magicancestry.statespace
    nfgl = length(statespace["haploindex"])
    d2g = rulediplo2geno(nfgl,ismtx=true)
    magicancestry.genoprob=[[round.(i*d2g,digits=posteriordigits) for i=prob]
        for prob=magicancestry.diploprob]
end

function prior_genoindex(nfgl::Integer)
    [[i,j] for i=1:nfgl for j=1:i]
end

function prior_diploindex(nfgl::Integer)
    [[i,j] for i=1:nfgl for j=1:nfgl]
end

function rulediplo2geno(nfgl::Integer;ismtx::Bool=true)    
    genotypes = prior_genoindex(nfgl)    
    diplotypes = prior_diploindex(nfgl)
    rulediplo2geno(diplotypes,genotypes; ismtx)    
end

function rulediplo2geno(inputdiplotypes::AbstractVector, genotypes::AbstractVector;ismtx::Bool=true)        
    issubset(genotypes,inputdiplotypes) || error("ancestral genotypes is not a subset of ancestral diplotypes")
    issubset(reverse.(genotypes),inputdiplotypes) || error("ancestral genotypes is not a subset of ancestral diplotypes")
    allunique(genotypes) || error("ancestral genotypes are not all unique")    
    diplotypes = [in(i, genotypes) ? i : reverse(i) for i = inputdiplotypes]
    genorule= Dict(genotypes .=> 1:size(genotypes,1))
    d2g=[get(genorule,i,0) for i=diplotypes]
    if ismtx
        d2gmtx= spzeros(length(diplotypes),length(genotypes))
        for i=eachindex(d2g)
            d2gmtx[i,d2g[i]]=1
        end
        d2gmtx
    else
        d2g
    end
end

"""
    sethaploprob!(magicancestry::MagicAncestry)

set magicancestry.haploprob from magicancestry.diploprob.

"""
function sethaploprob!(magicancestry::MagicAncestry; posteriordigits=4)
    if isnothing(magicancestry.diploprob)
        error("diploprob is required")
    end
    statespace=magicancestry.statespace
    d2h= rulediplo2haplo(statespace)
    magicancestry.haploprob=[[round.(i*d2h,digits=posteriordigits) for i=prob]
        for prob=magicancestry.diploprob]
end

function rulediplo2haplo(statespace::AbstractDict)
    nfgl = length(statespace["haploindex"])
    diplotypes = prior_diploindex(nfgl)
    nrow =size(diplotypes,1)
    d2h= zeros(nrow,nfgl)
    for i=1:nrow for j=1:2
            d2h[i,diplotypes[i][j]] += 0.5
        end
    end
    sparse(d2h)
end
function setinbredcoef!(magicancestry::MagicAncestry; posteriordigits)
    noff  = size(magicancestry.magicped.offspringinfo,1)
    nfgl = length(magicancestry.statespace["haplotype"])
    ibdcols =  [(i-1)nfgl+i for i in 1:nfgl]
    nchr = length(magicancestry.markermap)
    if !isnothing(magicancestry.diploprob)        
        magicancestry.inbredcoef = [begin 
            nsnp = size(magicancestry.markermap[chr],1)
            coef = zeros(nsnp, noff)
            prob = magicancestry.diploprob[chr]
            for off in 1:noff
                coef[:,off] .= round.(Vector(sum(view(prob[off],:,ibdcols),dims=2)[:,1]), digits=posteriordigits)
            end
            coef
        end for chr in 1:nchr]
    elseif !isnothing(magicancestry.viterbipath)        
        magicancestry.inbredcoef = [[in(i,ibdcols) for i in magicancestry.viterbipath[chr]] for chr in 1:nchr]
    else
        error("diploprob or viterbipath is required")
    end
end

function rulegeno2haplo(statespace::AbstractDict)
    nfgl = length(statespace["haploindex"])
    genotypes = prior_genoindex(nfgl)
    nrow =size(genotypes,1)
    g2h= zeros(nrow,nfgl)
    for i=1:nrow for j=1:2
            g2h[i,genotypes[i][j]] += 0.5
        end
    end
    sparse(g2h)
end

"""
    savecondprob(b(outputfile::AbstractString, magicancestry::MagicAncestry;
        probtype="haploprob"),workdir=pwd())

save probability field of magicancestry into a csv file or magicancestry into a jld2 file.

The output is magicancestry.haploprob, magicancestry.genoprob, magicancestry.diploprob or magicancestry for
probtype = "haploprob", "genoprob", or "diploprob", respectively.

"""
function savecondprob(outputfile::AbstractString,magicancestry::MagicAncestry;
    probtype::AbstractString="haploprob",
    delim::AbstractChar=',', workdir::AbstractString=pwd())
    if !(probtype in ["haploprob", "genoprob", "diploprob"])
        msg = string("unknown probtype: ", probtype)
        msg  = msg*". probtype must be \"haploprob\",\"genoprob\", or \"diploprob\"."
        error(msg)
    end
    if probtype == "haploprob"
        condprob = magicancestry.haploprob
    elseif probtype == "genoprob"
        condprob = magicancestry.genoprob
    elseif probtype == "diploprob"
        condprob = magicancestry.diploprob
    else
        error("unknown saving probaiblity type: ", probtype)
    end
    if isnothing(condprob)
        @warn string(prob, " not exist")
        return -1
    end
    outputfile2 = getabsfile(workdir,outputfile)
    open(outputfile2, "w") do io
        initial = nothing
        savecondprob_ind(io, magicancestry; delim, initial,probtype)
    end
    outputfile
end


function ancestrycall(condprob::AbstractVector;minprob::Real= 0.5)
    [begin
        chrprob = condprob[ch]
        nsnp = size(chrprob[1],1)
        chrbest= [begin
            val,index0 = findmax(chrprob[j][i,:])
            val > minprob ? index0 : missing
        end for j in eachindex(chrprob),  i=1:nsnp]
        permutedims(chrbest)
    end for ch in eachindex(condprob)]
end

function calnumrecom(offancestry::AbstractVector; isperchrom::Bool=false)
    # offancestry[chr][m,off]: ancesry for off at marker m
    nrecom=[[length(splitindex(collect(skipmissing(i))))-1
        for i=eachcol(a)] for a=offancestry]        
    nrecom2 = reduce(hcat, nrecom)
    if isperchrom
        nrecom2
    else
        sum(nrecom2,dims=2)[:,1]        
    end
end

function calnumrecom(magicancestry::MagicAncestry;
    isperchrom::Bool=false,
    minprob::Real= 0.7)
    prob  = magicancestry.genoprob
    prob2 = isnothing(prob) ? magicancestry.haploprob : prob
    if isnothing(prob2)
        bestancestry = magicancestry.viterbipath
    else
        bestancestry= ancestrycall(prob2,minprob=minprob)
    end
    nrecom = calnumrecom(bestancestry; isperchrom)
    nrecom
end

function submagicancestry(magicancestry::MagicAncestry;
    chrsubset::Union{Nothing,AbstractRange,AbstractVector}=nothing)
    nchr=length(magicancestry.markermap)
    if isnothing(chrsubset)
        chrsubset2 = 1:nchr
    else
        chrsubset2 = chrsubset[chrsubset.<= nchr]
        issetequal(chrsubset2, chrsubset) || @warn string("chrsubset is reset to ",chrsubset2)
    end
    issetequal(chrsubset2, 1:nchr) && return deepcopy(magicancestry)
    magicped = deepcopy(magicancestry.magicped)
    markermap = magicancestry.markermap[chrsubset2]
    foundergeno = magicancestry.foundergeno[chrsubset2]
    statespace =  deepcopy(magicancestry.statespace)
    if isnothing(magicancestry.viterbipath)
        viterbipath = nothing
    else
        viterbipath = magicancestry.viterbipath[chrsubset2]
    end
    if isnothing(magicancestry.diploprob)
        diploprob = nothing
    else
        diploprob = magicancestry.diploprob[chrsubset2]
    end
    if isnothing(magicancestry.genoprob)
        genoprob = nothing
    else
        genoprob = magicancestry.genoprob[chrsubset2]
    end
    if isnothing(magicancestry.haploprob)
        haploprob = nothing
    else
        haploprob = magicancestry.haploprob[chrsubset2]
    end
    if isnothing(magicancestry.inbredcoef)
        inbredcoef = nothing
    else
        inbredcoef = magicancestry.inbredcoef[chrsubset2]
    end
    loglike = magicancestry.loglike[chrsubset2, :]
    misc = deepcopy(magicancestry.misc)
    MagicAncestry(magicped,markermap,foundergeno,statespace,
        viterbipath,diploprob,genoprob,haploprob,inbredcoef, loglike,misc)
end

function mergeancestry(ancestryls::Vector{MagicAncestry})
    isempty(ancestryls) && @error string("ancestryls is empty")
    length(ancestryls) == 1 && return only(ancestryls)
    magicped = deepcopy(ancestryls[1].magicped)
    markermap = reduce(vcat,[i.markermap for i in ancestryls])
    foundergeno = reduce(vcat,[i.foundergeno for i in ancestryls])
    statespace = deepcopy(ancestryls[1].statespace)
    viterbils = reduce(vcat,[i.viterbipath for i in ancestryls])    
    viterbils = viterbils[.!isnothing.(viterbils)]
    if length(viterbils)==0
        # isempty(viterils) results in false!
        viterbipath = nothing
    else
        viterbipath = viterbils
    end    
    probls = reduce(vcat,[i.diploprob for i in ancestryls])
    deleteat!(probls,isnothing.(probls))    
    diploprob = length(probls)==0 ? nothing : probls
    probls = reduce(vcat,[i.genoprob for i in ancestryls])
    deleteat!(probls,isnothing.(probls))    
    genoprob = length(probls)==0 ? nothing : probls
    probls = reduce(vcat,[i.haploprob for i in ancestryls])
    deleteat!(probls,isnothing.(probls))    
    haploprob = length(probls)==0 ? nothing : probls
    coefls = reduce(vcat,[i.inbredcoef for i in ancestryls])    
    deleteat!(coefls,isnothing.(coefls))    
    inbredcoef = length(coefls)==0 ? nothing : coefls
    loglike = reduce(vcat,[i.loglike for i in ancestryls])
    misc = deepcopy(ancestryls[1].misc)
    MagicAncestry(magicped,markermap,foundergeno,statespace,
        viterbipath,diploprob,genoprob,haploprob,inbredcoef,loglike,misc)
end

function deloutlier!(magicancestry::MagicAncestry)
    isoutlier = magicancestry.magicped.offspringinfo[!,:isoutlier]
    if any(isoutlier)
        @info string("delete ",sum(isoutlier), " outlier offspring: ", magicancestry.magicped.offspringinfo[isoutlier, :individual])
    else
        @info string("No outlier offspring")
    end
    isnonoutlier = .!isoutlier
    magicancestry.magicped.offspringinfo = magicancestry.magicped.offspringinfo[isnonoutlier, :]
    magicancestry.loglike = magicancestry.loglike[:,isnonoutlier]
    for chr in eachindex(magicancestry.markermap)
        if !isnothing(magicancestry.diploprob)
            magicancestry.diploprob[chr] = magicancestry.diploprob[chr][isnonoutlier]
        end
        if !isnothing(magicancestry.genoprob)
            magicancestry.genoprob[chr] = magicancestry.genoprob[chr][isnonoutlier]
        end
        if !isnothing(magicancestry.haploprob)
            magicancestry.haploprob[chr] = magicancestry.haploprob[chr][isnonoutlier]
        end    
        if !isnothing(magicancestry.inbredcoef)
            magicancestry.inbredcoef[chr] = magicancestry.inbredcoef[chr][:, isnonoutlier]
        end    
        if !isnothing(magicancestry.viterbipath)
            magicancestry.viterbipath[chr] = magicancestry.viterbipath[chr][:, isnonoutlier]
        end
    end
    magicancestry
end


function thinmagicancestry!(magicancestry::MagicAncestry;thincm::Real=0.0,isdeloutlier::Bool=true)
    if isdeloutlier
        deloutlier!(magicancestry)
    end
    if thincm < 0         
        return magicancestry
    end
    for chr in eachindex(magicancestry.markermap)
        # always keep first and last marekrs
        d = diff(magicancestry.markermap[chr][!,:poscm])
        pushfirst!(d,100)
        d[end] = 100
        snps = findall(d .> thincm)
        magicancestry.markermap[chr] = magicancestry.markermap[chr][snps,:]
        magicancestry.foundergeno[chr] = magicancestry.foundergeno[chr][snps,:]
        if !isnothing(magicancestry.haploprob)
            for off in eachindex(magicancestry.haploprob[chr])
                magicancestry.haploprob[chr][off] =  magicancestry.haploprob[chr][off][snps,:]
            end
        end
        if !isnothing(magicancestry.genoprob)
            for off in eachindex(magicancestry.genoprob[chr])
                magicancestry.genoprob[chr][off] =  magicancestry.genoprob[chr][off][snps,:]
            end
        end
        if !isnothing(magicancestry.diploprob)
            for off in eachindex(magicancestry.diploprob[chr])
                magicancestry.diploprob[chr][off] =  magicancestry.diploprob[chr][off][snps,:]
            end
        end
        if !isnothing(magicancestry.viterbipath)    
            magicancestry.viterbipath[chr] = magicancestry.viterbipath[chr][snps,:]
        end
        if !isnothing(magicancestry.inbredcoef)
            magicancestry.inbredcoef[chr] = magicancestry.inbredcoef[chr][snps,:]
        end
    end    
    magicancestry
end


function thinmagicancestry(ancestryfile::AbstractString;
    thincm::Real=0, 
    isdeloutlier::Bool=true, 
    workdir::AbstractString=pwd(),
    outext::Union{Nothing,AbstractString}=".csv.gz",
    outstem::AbstractString=first(MagicBase.split_allext(basename(ancestryfile))),
    verbose::Bool=true)
    magicancestry = readmagicancestry(ancestryfile;workdir)   
    if verbose
        nsnps = sum(size.(magicancestry.markermap,1))
        @info string("#snps = ",nsnps, " for inputfile =",ancestryfile)
    end
    MagicBase.thinmagicancestry!(magicancestry; thincm,isdeloutlier)
    if isnothing(outext)
        outext = last(split_allext(ancestryfile))
    end
    outfile = string(outstem,"_thincm",string(thincm),outext)
    outfile2 = getabsfile(workdir,outfile)
    if verbose
        nsnps = sum(size.(magicancestry.markermap,1))
        @info string("#snps = ",nsnps, " for outputfile =",outfile)
    end
    savemagicancestry(outfile2,magicancestry)
end