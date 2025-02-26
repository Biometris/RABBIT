
# https://discourse.julialang.org/t/how-to-get-current-julia-process-memory-usage/41734/4
function memoryuse()
    if Sys.isunix()
        open("/proc/$(getpid())/stat") do io
            s = read(io, String)
            tryparse( Int, split(s)[23])
        end
    elseif Sys.iswindows()
        s = String(read(`tasklist /FI "PID eq $(getpid())" /NH /FO csv`))
        s = replace(s, "\"" => "")
        # @info s          
        m = join(strip.(split(s, ","))[5:end])
        m = replace(m, "."=>"",","=>"", "K" => "1000", "M" => "1000000", "G" => "1000000000")        
        try 
            tryparse.(Float64, split(m)) |> prod |> Int
        catch err
            @warn err
            @warn "Could not get mememoryuse."
            -1
        end
    else 
        -1
    end    
end

function pkgversion_st(pkgname::AbstractString)    
    io = IOBuffer()
    Pkg.status(pkgname;io)
    st = String(take!(io))
    st2 = split(st," ")
    pos = findfirst(st2 .== pkgname)
    if isnothing(pos)
        @warn string(pkgname, " does not exist or are not installed")
        nothing
    else
        VersionNumber(strip(replace(st2[pos+1],"v"=>"")))
    end
end

function pkgversion_st(m::Module)
    pkgname = string(m)
    pkgversion_st(pkgname)    
end

function pkgversion_toml(m::Module)
    repodir = pkgdir(m)
    tomlfile = joinpath(repodir, "Project.toml")
    isfile(tomlfile) || @error string(tomlfile, "does not exist") 
    VersionNumber(Pkg.TOML.parsefile(tomlfile)["version"])
end

function printconsole(io::Union{Nothing,IO},verbose::Bool,msg::AbstractString)
    verbose && @info(msg)
    if !isnothing(io)
        write(io,string(msg,"\n"))
        flush(io)
    end
end


# PkgVersion
# https://discourse.julialang.org/t/how-to-find-out-the-version-of-a-package-from-its-module/37755/5
function printpkgst(io::Union{Nothing,IO},verbose::Bool,pkgname::AbstractString)
    verbose && Pkg.status(pkgname)    
    if !isnothing(io)
        Pkg.status(pkgname;io)
    end
end

function getabsfile(workdir::AbstractString,file::AbstractString)
    isdir(workdir) || error(string("workdir=", workdir, " is not a directory"))
    a = splitpath(file)
    if length(a) > 1
        file
    else
        joinpath(workdir,file)
    end
end

function splitindex(f::Function,A::AbstractVector)
    size(A,1)==1 && return [1:1]
    res=Vector{typeof(1:1)}()
    i0=1
    for i=2:size(A,1)
        if !f(A[i-1],A[i])
            push!(res,i0:i-1)
            i0=i
        end
    end
    push!(res,i0:size(A,1))
    res
end

function splitindex(A::AbstractVector)
    f(x,y)= if ismissing(x) 
        ismissing(y)
    else
        if ismissing(y)
            false
        else
            x==y
        end
    end        
    splitindex(f,A)
end

function logsumexp(a::AbstractVector)
    isempty(a) && return 0.0
    off, ind = findmax(a)
    b = exp.(a .- off)
    b[ind] -= 1.0
    # using KahanSummation
    # log1p(sum_kbn(b))+off
    log1p(sum(b))+off
end

function readdataframe(filename::AbstractString;
    delim::AbstractChar= ',',
    commentstring::AbstractString="##",    
    missingstring=["NA","missing"],
    workdir::AbstractString=pwd())
    ext = last(splitext(filename))
    myopen = ext == ".gz" ? GZip.open : open
    myopen(getabsfile(workdir,filename), "r") do io
        lines = readlines(io,keep=true)
        b = startswith.(lines,commentstring)
        deleteat!(lines,b)        
        length(lines) >= 1 || @error "empty file"        
        df = CSV.read(IOBuffer(join(lines)),DataFrame;delim, missingstring)                
        df
    end
end

function readmultitable(filename::AbstractString;
    delim::AbstractChar=',',
    commentstring::AbstractString="##",
    missingstring=["NA","missing"], # not adopted when isparse is false
    isparse::Bool= true, 
    isordereddict::Bool = true)
    ext = last(splitext(filename))
    myopen = ext == ".gz" ? GZip.open : open
    myopen(filename, "r") do io
        lines = readlines(io,keep=isparse)        
        b = startswith.(lines,commentstring)
        deleteat!(lines,b)
        length(lines) >= 1 || @error "empty file"
        regionkey = first(split(lines[1],delim))
        seg = findall(startswith.(lines,regionkey))
        push!(seg,length(lines)+1)
        mydict = isordereddict ? OrderedDict : Dict    
        mydict([begin
            i1=seg[i]
            i2 =seg[i+1]
            tabtitle = split(lines[i1],delim)            
            length(tabtitle) >= 2 || error("tabe name not exist")
            tabkey = strip(chomp(tabtitle[2]),'\"')
            tabval = view(lines, i1+1:i2-1)                        
            if isparse
                # bug in CSV.read, a cell with value of  a very very long string will be truncated!!
                # for example, sparse probability matrix in magicancestry, 
                val = CSV.read(IOBuffer(join(tabval)),DataFrame;missingstring)                
            else
                strval = reduce(hcat, split.(chomp.(tabval),delim))          
                mtx =  permutedims(strip.(strval[:,2:end],'\"'))
                cols = strip.(strval[:,1],'\"')
                val = DataFrame(mtx,cols)
            end
            tabval .= ""                        
            tabkey => val
        end for i=1:length(seg)-1])
    end
end


function set_logfile_begin(logfile::Union{Nothing,AbstractString,IO},
    workdir::AbstractString,
    funname::AbstractString;
    verbose::Bool=true,
    delim::Union{Nothing,AbstractString}="-")
    if isnothing(logfile)
        io=nothing
    else
        if isa(logfile, AbstractString)
            io=open(MagicBase.getabsfile(workdir,logfile), "w")
            msg = string(funname, ", logfile=", logfile,
                ", ", round(Dates.now(),Dates.Second))
            printconsole(io,verbose,msg)
        else
            # isa(io, IO)
            io=logfile
        end
    end
    if !isnothing(delim)
        nleft = (76 - length(funname)) ÷ 2
        nright = 76 - nleft -length(funname)
        MagicBase.printconsole(io,verbose,string(repeat(delim,nleft),
            funname,repeat(delim,nright)))
    end
    io
end

function set_logfile_end(logfile::Union{Nothing,AbstractString,IO},
    io::Union{Nothing,AbstractString,IO},
    tused::Real,
    funname::AbstractString;
    verbose::Bool=true,
    delim::Union{Nothing,AbstractString}="-")
    msg= string("End, ", round(Dates.now(),Dates.Second),", tused = ",
        tused, " seconds by ",funname)
    MagicBase.printconsole(io,verbose,msg)
    isnothing(delim) || MagicBase.printconsole(io,verbose,repeat(delim,76))
    isa(logfile, AbstractString) && close(io)
end

function split_allext(filename::AbstractString)
    list_rabbitext = [".vcf",".csv",".tar",".gz",".png",".log",".jl",".txt",".tsv",".xlsx",".pdf",".zip",".7z"]
    filebase,ext = splitext(filename)
    isempty(ext) && return (filebase,ext)
    while true    
        filebase,ext2 = splitext(filebase)
        if isempty(ext2)  
            break
        else
            if ext2 in list_rabbitext
                ext = ext2 * ext
            else
                filebase *=ext2
                break
            end
        end        
    end
    filebase,ext
end

function create_tar(dir::AbstractString,outfile_tar::AbstractString)
    if !occursin(r".tar$",outfile_tar)
        @warn string(outfile_tar, " is not .tar")
    end
    isdir(dir) || @error string(dir, " is not a directory")
    Tar.create(dir, outfile_tar)
end

function create_targz(dir::AbstractString,outfile_targz::AbstractString)
    if !occursin(r".tar.gz$",outfile_targz)
        @warn string(outfile_targz, " is not .tar.gz")
    end
    isdir(dir) || @error string(dir, " is not a directory")
    open(outfile_targz, "w") do tar_gz    
        tar = GzipCompressorStream(tar_gz)
        try
            Tar.create(dir, tar)
        finally
            close(tar)
        end    
    end
    outfile_targz
end

function extract_tar(tarfile::AbstractString, outdir::AbstractString)
    if !occursin(r".tar$",tarfile)
        @warn string(tarfile, " is not .tar")
    end
    isdir(outdir) || @error string(outdir, " is not a directory")
    Tar.extract(tarfile,outdir)
end

function extract_targz(targzfile::AbstractString, outdir::AbstractString)
    if !occursin(r".tar.gz$",targzfile)
        @warn string(targzfile, " is not .tar.gz")
    end
    isdir(outdir) || @error string(outdir, " is not a directory")
    tar_gz = open(targzfile)
    try
        tar = GzipDecompressorStream(tar_gz)
        try
            Tar.extract(tar,outdir)
        finally
            close(tar)
        end
    finally
        close(tar_gz)
    end
end

function read_headlines(file::AbstractString, nline::Integer)
    nline <= 0 && return ""
    isgz = last(splitext(file)) == ".gz"
    myopen = isgz ? GZip.open : open
    myopen(file,"r") do io
        lines = ""
        for i in 1:nline
            eof(io) && break
            line = readline(io,keep=true)
            lines *= line
        end
        rstrip(lines,['\n'])
    end
end

function findlastcomment(file::AbstractString; commentstring::AbstractString="##")
    isgz = last(splitext(file)) == ".gz"
    myopen = isgz ? GZip.open : open
    myopen(file,"r") do io
       commentline = 0
       while !eof(io)
           word = readline(io)
           if isempty(word) 
                break
           else
                startswith(word,commentstring) ? (commentline+=1) : break
           end
       end
       commentline
    end
end


function info_argnames(@nospecialize(f))
    mls = methods(f)
    for i in eachindex(mls)
        posargls = Base.method_argnames(mls[i])
        popfirst!(posargls)
        keyargls = Base.kwarg_decl(mls[i])
        if length(keyargls) <= 5
            msg = string(f, " method ",i, ": posarg=",posargls,", keyarg=",keyargls)
        else
            msg = string(f, " method ",i, "\n posarg=",posargls,"\n keyarg=",keyargls)
        end
        @info msg
    end
    nothing
end

function metroplis1d(loglfun::Function;
    xstart::Real=0.0, temperature::Real=1.0,
    constraint::Function = x->true,
    stepsize::Real=1.0, nstep::Integer=1)
    xnow = xstart
    logl = loglfun(xnow)
    chain = [[xnow, logl,NaN]]
    for it in 1:nstep
        xprop = rand(Normal(xnow, stepsize))
        if constraint(xprop)
            proplogl = loglfun(xprop)
            if temperature < 10^(-5.0)
                cond = proplogl > logl
            else
                cond = rand() < min(1.0, exp((proplogl-logl)/temperature))
            end
            if cond
                xnow, logl = xprop, proplogl
            end
        else
            cond = false
        end
        push!(chain, [xnow, logl, cond])
    end
    chain
end