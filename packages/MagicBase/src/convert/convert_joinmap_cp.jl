

# Generate jointmap and rabbit input from jointmap body geno datafile.
# TSV inputfile: a dataframe file with columns: marker, seggregation, offspring1, offspring2, ...
# joionmap2rabbit works only for CP population 

function joionmap2rabbit(cpfile::AbstractString;
    delcols::Union{Nothing,AbstractVector}=nothing, 
    delrows::Union{Nothing,AbstractVector}=nothing, 
    parentid::AbstractVector = ["P1","P2"],        
    familyid::AbstractString = "CP",
    outstem = "outstem",    
    logfile::Union{AbstractString,IO}= string(outstem,".log"),
    workdir::AbstractString = pwd(),
    commentstring::AbstractString="##",
    verbose::Bool=true)
    starttime = time()
    io = MagicBase.set_logfile_begin(logfile, workdir, "rabbit_cp_input"; verbose)
    isa(logfile, AbstractString) && printpkgst(io,verbose,"MagicBase")
    msg = string("list of args: \n",
        "cpfile = ", cpfile, "\n",        
        "delcols = ", delcols, "\n",        
        "delrows = ", delrows, "\n",      
        "parentid = ", parentid, "\n",   
        "familyid = ", familyid, "\n",   
        "outstem = ", outstem, "\n",       
        "logfile = ", logfile, "\n",       
        "workdir = ", workdir, "\n",
        "verbose = ", verbose)
    printconsole(io,verbose,msg)
    @warn "current version works only CP population"
    cpgeno = MagicBase.read_joinmap_cpfile(cpfile; delcols,delrows,io,workdir);
    genofile,pedfile = rabbit_cp_input!(cpgeno; parentid,familyid, outstem,io,workdir,commentstring)
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, io, tused,"rabbit_cp_input"; verbose)
    genofile,pedfile 
end

function read_joinmap_cpfile(cpfile::AbstractString; 
    delcols::Union{Nothing,AbstractVector}=nothing, 
    delrows::Union{Nothing,AbstractVector}=nothing, 
    io::Union{Nothing,IO}=nothing,
    workdir::AbstractString=pwd(),    
    verbose::Bool=true)    
    cpgeno = CSV.read(getabsfile(workdir,cpfile),DataFrame)    
    printconsole(io,verbose, "beginning input geno: ")
    eg_geno = cpgeno[1:min(2,size(cpgeno,1)),1:min(8,size(cpgeno,2))]    
    verbose && @info eg_geno
    isnothing(io) || CSV.write(io,eg_geno;append=true,writeheader=true)
    if !isnothing(delcols) 
        msg = string("delete columns: ",join(names(cpgeno,delcols),","))
        printconsole(io,verbose,msg)
        select!(cpgeno,Not(delcols))
    end
    isnothing(delrows) || deleteat!(cpgeno,delrows)
    col12=names(cpgeno,1:2)
    newcol12 = ["marker","segregation"]
    msg = string("rename first 2 columns ",col12, " => ",newcol12)
    printconsole(io,verbose,msg)
    rename!(cpgeno,Dict(col12 .=> newcol12))
    # check segregation types
    segrls = ["<abxcd>","<efxeg>","<hkxhk>","<lmxll>","<nnxnp>"]
    d = setdiff(unique(cpgeno[!,:segregation]),segrls)    
    if !isempty(d) 
        msg = string("seggreation types  ",join(d[1:min(10,length(d))],","),", not ",join(segrls,",")," for CP")
        printconsole(io,false,msg)        
        @error msg
    end    
    # check genotypes of a few beginning individuals
    indls = names(cpgeno,3:min(7,size(cpgeno,2)))    
    # missing allele codes: ".","u", "-"
    alleleset = [".","u", "-", "a","b","c","d","e","f","g","h","k","l","m","n","p"]    
    for col in indls
        als = lowercase.(unique(split(join(unique(cpgeno[!,col])),"")))
        adiff = setdiff(als,alleleset)
        if !isempty(adiff)
            msg = string("unknown alleles in col=",col,": ", adiff[1:min(10,length(adiff))])
            printconsole(io,false,"ERROR: "*msg)        
            @error msg
        end
    end            
    printconsole(io,verbose, "beginning transformed geno: ")
    eg_geno = cpgeno[1:min(2,size(cpgeno,1)),1:min(8,size(cpgeno,2))]    
    verbose && @info eg_geno
    isnothing(io) || CSV.write(io,eg_geno;append=true,writeheader=true)
    msg = string("#marker=",size(cpgeno,1),", #offspring=",size(cpgeno,2)-2)
    msg *= string(", 1st-offspring=",first(indls))    
    printconsole(io,verbose,msg)
    cpgeno
end

function rabbit_cp_input!(cpgeno::DataFrame;
    parentid::AbstractVector = ["P1","P2"],        
    familyid::AbstractString = "CP",
    outstem = "outstem",    
    io::Union{Nothing,IO}=nothing,
    workdir::AbstractString = pwd(),
    commentstring::AbstractString="##",
    verbose::Bool=true)
    # insert parent genotypes based segregation
    p1p2 = [string.(split(replace(i,"<"=>"",">"=>""),"x")) for i in cpgeno[!,:segregation]]
    insertcols!(cpgeno,3,Symbol(parentid[1]) => first.(p1p2),Symbol(parentid[2]) => last.(p1p2))
    # convert alleles to 0, 1, . (missing)
    # missing allele codes: ".","u", "-"
    for row in eachrow(cpgeno)
        segr = row[:segregation]
        allelerule = ["-"=>".","u"=>"."]
        alleles = sort(unique(split(join(split(replace(segr,"<"=>"",">"=>""),"x")),"")))
        append!(allelerule,(alleles .=> string.(0:(length(alleles)-1))))
        row[3:end] .= [replace(lowercase(i),allelerule...) for i in row[3:end]]
    end
    reshead = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL","FILTER", "INFO", "FORMAT"]
    append!(reshead,names(cpgeno,3:size(cpgeno,2)))
    resleft = Matrix{String}(undef, size(cpgeno,1), 9)
    resleft .= "."
    # VCFv4.3 requires that REF must be one of A, T, G, C, N
    resleft[:,4] .= "N" 
    resleft[:,3] .= cpgeno[!,:marker]
    # resleft[:,8] .= ["SEGR="*i for i in cpgeno[!,:segregation]]
    resleft[:,9] .= "GT"
    res = Matrix{String}(cpgeno[!,3:end])    
    empty!(cpgeno)
    gset = unique(res)
    grule =gset .=> join.(split.(gset,""),"/")    
    replace!(res,grule...)
    res = hcat(resleft,res)
    # print beginning vcf
    nind = min(4,length(reshead)-9)
    cols  = vcat(1:3,8:9,10:9+nind)
    printconsole(io,verbose, string("beginning vcf geno of cols= ",cols))    
    eg_geno = DataFrame(res[1:2,cols],reshead[cols])
    verbose && @info eg_geno
    isnothing(io) || CSV.write(io,eg_geno;append=true,writeheader=true)
    # save results
    outfile = outstem*"_geno.vcf.gz"
    GZip.open(getabsfile(workdir,outfile),"w") do outio
        cc = commentstring
        write(outio, cc*"fileformat=VCFv4.3\n")
        filedate = string(cc*"##filedate=",string(now()))
        write(outio, filedate,"\n")
        write(outio, cc*"source=RABBIT\n")        
        msg = cc*"FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        write(outio, msg,"\n")        
        msg = cc*"INFO=<ID=SEGR,Number=1,Type=String,Description=\"jointmap segregation type for cp population\">"
        write(outio, msg,"\n")
        write(outio, join(reshead,"\t"),"\n")
        writedlm(outio,res,'\t')
    end
    offls = reshead[12:end]
    msg = string("#marker=",size(res,1),", #parent=2, #offspring=",length(offls))     
    printconsole(io,verbose,msg)
    printconsole(io,verbose,string("output genofile: ",outfile))
    # save pedfile
    ped = MagicBase.parsedesign("2ril-self0";popid = familyid)
    magicped = formmagicped(ped,length(offls); isfglexch=false)
    MagicBase.setfounderid!(magicped,parentid)
    MagicBase.setoffspringid!(magicped,offls)
    pedfile = outstem*"_ped.csv"
    savemagicped(getabsfile(workdir,pedfile),magicped)
    printconsole(io,verbose,string("output pedfile: ",pedfile))
    outfile,pedfile
end
