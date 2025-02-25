function check_infiles(genofile::AbstractString,pedinfo::Union{Integer,AbstractString}; 
    isbreedped::Union{Nothing,Bool}=nothing,
    ismagicsimulate::Bool=false,
    commentstring = "##",
    workdir::AbstractString=pwd(),
    io::Union{Nothing,IO}=nothing,
    verbose::Bool=true)
    # check genofile extension and exist
    fileext = last(MagicBase.split_allext(genofile))
    if in(fileext,[".vcf",".vcf.gz"])         
        msg = string("check genofile extension = ",fileext, ": OK")            
        printconsole(io,verbose,msg)        
    elseif in(fileext,[".csv",".csv.gz"]) 
        msg = string("check genofile extension = ",fileext, ": CSV-format genofile is not recommended")
        printconsole(io,false,"Warning: "*msg)
        @warn msg        
    else
        msg = string("check genofile extension = ",fileext, ": FAILED. Must be .vcf or .vcf.gz")
        printconsole(io,false,"Error: "*msg)
        error(msg)
    end
    genofile2 = getabsfile(workdir,genofile)
    if isfile(genofile2) 
        msg =  string("check genofile existence: YES")
        printconsole(io,verbose,msg)        
    else
        msg = string(genofile2, " does not exist")
        printconsole(io,false,"Error: "*msg)
        error(msg)
    end
    # check pedinfo and exist
    if isa(pedinfo,Integer)
        if pedinfo < 0 
            msg = string("Integer pedinfo = ",pedinfo, ": FAILED. Must be non-negative for the number of founders")
            printconsole(io,false,"Error: "*msg)
            error(msg)
        else
            msg = string("Integer pedinfo = ",pedinfo, ": OK")
            printconsole(io,verbose,msg)                    
        end        
    elseif isa(pedinfo,AbstractString)
        # pedinfo = "nfounder=19||ibd=0.9852||j1122=5.469||j1211=0.0919||j1213=0.0537||j1222=0.0919||j1232=0.0537"
        # if pedinfo  is juncdist designinfo, it may return fileext such as .0537
        fileext = last(MagicBase.split_allext(pedinfo))
        if isempty(fileext) || occursin("=",pedinfo)  
            # designcode
            designinfo = parsedesign(pedinfo) 
            if !isnothing(designinfo) 
               msg = string("check pedinfo as stringcode: OK")
               printconsole(io,verbose,msg)        
            end            
        else
            # pedigree file
            if in(fileext,[".csv",".csv.gz"]) 
                msg = string("check pedfile extension = ",fileext, ": OK")            
                printconsole(io,verbose,msg)        
            else
                msg = string("check pedfile extension = ",fileext, ": FAILED. Must be .csv or .csv.gz")
                printconsole(io,false,"Error: "*msg)
                error(msg)
            end
            pedfile2 = getabsfile(workdir,pedinfo)
            if isfile(pedfile2) 
                msg = string("check pedfile existence: YES")
                printconsole(io,verbose,msg)        
            else
                printconsole(io,false,"Error: "*msg)
                error(string(pedfile2, " does not exist"))
            end            
        end
    else
        msg = string("wrong pedinfo = ",pedinfo)
        printconsole(io,false,"Error: "*msg)
        error(msg)
    end
    # check consistency betwen sampleid of genofile and pedinfo
    if in(last(MagicBase.split_allext(genofile)),[".csv",".csv.gz"]) 
        samples = csv_get_samples(genofile; commentstring, workdir)    
    else
        samples = vcf_get_samples(genofile; commentstring, workdir)    
    end
    check_allunique(samples,"samples";io, verbose)    
    if isa(pedinfo,Integer)
        nfounder = pedinfo 
        if length(samples) >= nfounder
            msg = string("check founders: OK")            
            printconsole(io,verbose,msg)
        else
            msg = string("check founders: FAILED. #samples = ",length(samples), " in genofile < #founders =",nfounder, " in pedinfo = ",pedinfo)
            printconsole(io,false,"Error: "*msg)
            error(msg)
        end
    elseif isa(pedinfo,AbstractString)
        fileext = last(MagicBase.split_allext(pedinfo))
        if isempty(fileext) || occursin("=",pedinfo)  
            designinfo = parsedesign(pedinfo)            
            if in(designinfo.designtype,[:commoncross,:juncdist])
                nfounder = getnfounder(designinfo)    
                if length(samples) >= nfounder
                    msg = string("check founders: OK if the intial ", nfounder, " samples = ", samples[1:nfounder], 
                        " are assumed to be founders for pedinfo = ",pedinfo)
                    printconsole(io,verbose,msg)        
                    # msg = string( "to confirm that founders = ", samples[1:nfounder])
                    # @warn msg
                    # printconsole(io,false,"Warning: "*msg)        
                    msg = string("check offspring: OK if the rest ", length(samples)-nfounder, " samples are assumed to be offspring")
                    printconsole(io,verbose,msg)        
                else
                    msg = string("check founders: FAILED. #samples = ",length(samples), " in genofile < #founders =",nfounder, " in pedinfo = ",pedinfo)
                    printconsole(io,false,"Error: "*msg)
                    error(msg)
                end                
            elseif designinfo.designtype == :breedcross
                founders = designinfo.founders
                check_allunique(founders,"founders";io, verbose)    
                d = setdiff(founders,samples)
                if isempty(d)                     
                    msg = string("check founders: OK.  All founders = ", founders, " specified by pedinfo =", pedinfo, " are also in samples in genofile")
                    printconsole(io,verbose,msg)        
                    msg = string("check offspring: OK if the rest ", length(samples)-length(founders), " samples are assumed to be offspring")
                    printconsole(io,verbose,msg)        
                else
                    msg = string("check founders: FAILED. Among ", length(founders), " founders, ",length(d), " founders = ", d,  " specified by pedinfo =", pedinfo, " are NOT in samples in genofile")
                    printconsole(io,false,"Error: "*msg)
                    error(msg)                    
                end
            else
                @error string("unknown designtype=",designtype, " for pedinfo = ",pedinfo)
            end    
        else            
            if isnothing(isbreedped)
                pedfile2 = getabsfile(workdir,pedinfo)   
                israbbitped = get_israbbitped(pedfile2; io)
                isbreedped = !israbbitped                         
                msg = string("reset isbreedped = ",isbreedped, ", pedfile=",pedinfo, " detected to be in ", isbreedped ? "breedped" : "RABBITped", " format")
                printconsole(io,verbose,msg)
            end
            if isbreedped                
                magicped = parsebreedped(pedinfo; outfile=nothing, commentstring, workdir)
            else
                magicped = readmagicped(pedinfo; commentstring, workdir)
            end
            founders = magicped.founderinfo[!,:individual]
            check_allunique(founders,"founders";io, verbose)    
            d = setdiff(founders,samples)            
            if isempty(d)                     
                msg = string("check founders: OK. All ", length(founders), " founders = ", founders, " in pedfile are also in samples in genofile")
                printconsole(io,verbose,msg)        
            else
                if length(d) == length(founders)                    
                    msg = string("check founders: FAILED. All ", length(d), " founders = ", d, " in pedfile but not in samples in genofile")
                    printconsole(io,false,"Error: "*msg)
                    error(msg)
                else
                    msg = string("check founders: PASS. Among ", length(founders), " founders, ", length(d), " founders = ", d, " in pedfile but not in samples in genofile. Ignore these founderrs and their offspring.")                    
                    printconsole(io,false,"Warning: "*msg)
                    @warn msg
                end
            end
            offspring = magicped.offspringinfo[!,:individual]
            check_allunique(offspring,"offspring";io, verbose)    
            if ismagicsimulate
                0
            else
                d = setdiff(offspring,samples)                        
                if isempty(d)                     
                    msg = string("check offspring: OK. All ", length(offspring), " offspring in pedfile are also in samples in genofile")
                    printconsole(io,verbose,msg)        
                else
                    msg = string("check offspring: PASS. Among ", length(offspring), " offspring, ignore ", length(d), " offspring = ", d, " in pedfile but not in samples in genofile")
                    printconsole(io,false,"Warning: "*msg)
                    @warn msg
                end
            end
            d = setdiff(samples, founders, offspring)
            if isempty(d)                     
                msg = string("check samples: OK. All ", length(samples), " samples in genofile are also in pedfile")
                printconsole(io,verbose,msg)        
            else
                msg = string("check samples: PASS. Among ", length(samples), " samples, ignore ", length(d), " samples in genofile but not in pedfile")
                printconsole(io,false,"Warning: "*msg)
                @warn msg
            end
        end
    else
        msg = string("wrong pedinfo = ",pedinfo)
        printconsole(io,false,"Error: "*msg)
        error(msg)
    end
    msg = string("Pass all checks for \n genofile = ",genofile,"\n pedinfo = ",pedinfo)
    if length(splitpath(genofile)) == 1 && length(splitpath(string(pedinfo))) == 1         
        msg *= string("\n workdir = ",workdir)
    end    
    printconsole(io,verbose,msg)        
end

function check_allunique(samples::AbstractVector, varname::AbstractString;
    io::Union{Nothing,IO}=nothing,
    verbose::Bool=true)
    if allunique(samples)        
        msg = string("check ", varname, " allunique: YES")            
        printconsole(io,verbose,msg)
    else
        ls = unique(samples)
        isdupe = [count(==(i)) > 1 for i in ls]
        msg = string("check ", varname, " allunique: FAILED. Duplicated ", varname, " = ",ls[isdupe])            
        printconsole(io,verbose,msg)
    end
end

function get_israbbitped(pedfile::AbstractString;io::Union{Nothing,IO}=nothing)
    lines = readlines(pedfile)
    b1 = occursin.(r"RABBIT, *designinfo",lines)
    b2 = occursin.(r"RABBIT, *offspringinfo",lines)
    israbbitped = any(b1) || any(b2)
    if israbbitped 
        if sum(b1) == 0 
            msg = string("designinfo is missing in RABBIT pedfile=",pedfile)
            printconsole(io,false,"Warning: "*msg)
            @warn msg
        elseif sum(b1) > 1
            msg = string("too many designinfo in RABBIT pedfile=",pedfile)
            printconsole(io,false,"Warning: "*msg)
            @warn msg
        end
        if sum(b2) == 0 
            msg = string("offspringinfo is missing in RABBIT pedfile=",pedfile)
            printconsole(io,false,"Warning: "*msg)
            @warn msg
        elseif sum(b2) > 1
            msg = string("too many offspringinfo in RABBIT pedfile=",pedfile)
            printconsole(io,false,"Warning: "*msg)
            @warn msg
        end
    end
    israbbitped
end

