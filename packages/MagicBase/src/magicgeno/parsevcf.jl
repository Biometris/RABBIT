
function parsevcf(genofile::AbstractString,pedinfo::Union{Integer,AbstractString};
    isfounderinbred::Bool=true,
    formatpriority::AbstractVector=["AD","GT"],
    isdelmultiallelic::Bool=true,    
    commentstring::AbstractString="##",
    outstem::Union{IO,AbstractString} = "outstem",
    outext::AbstractString = ".csv.gz",
    logfile::Union{Nothing, AbstractString,IO} = (isa(outstem, AbstractString) ? outstem : first(split_allext(genofile)))*"_parsevcf.log",
    workdir::AbstractString=pwd(),
    verbose::Bool=true)
    starttime = time()
    logio = MagicBase.set_logfile_begin(logfile, workdir, "parsevcf"; verbose,delim="-")
    inext = last(split_allext(genofile))
    inext in [".vcf",".vcf.gz"] || @error string("genofile ext must be .vcf or .vcf.gz")
    in_open = inext in [".vcf.gz"] ? GZip.open : open
    msg = string("list of args: \n",
        "genofile = ", genofile, "\n",
        "pedinfo = ", pedinfo, "\n",
        "isfounderinbred = ", isfounderinbred, "\n",
        "formatpriority = ", formatpriority, "\n",
        "isdelmultiallelic = ", isdelmultiallelic, "\n",
        "commentstring = ", commentstring, "\n",
        "outstem = ", outstem, "\n",
        "outext = ", outext, "\n",
        "logfile = ", logfile, "\n",
        "workdir = ", workdir, "\n",
        "verbose = ", verbose)
    printconsole(logio,verbose,msg)
    genofile2=getabsfile(workdir,genofile)
    lastcomment = findlastcomment(genofile2; commentstring)
    nheader = lastcomment + 1
    outext in [".csv",".csv.gz",".vcf",".vcf.gz"] || @error string("outfile ext must be in [.csv,.csv.gz,.vcf,.vcf.gz]")
    keepvcf = outext in [".vcf",".vcf.gz"]
    in_open(genofile2,"r") do inio
        if isa(outstem, AbstractString)
            out_open = outext in [".vcf.gz",".csv.gz"] ? GZip.open : open
            outfile = getabsfile(workdir,outstem*outext)
            out_open(outfile, "w") do outio
                parsevcf_io(inio, outio, logio, pedinfo,isfounderinbred,
                    formatpriority,isdelmultiallelic, commentstring, nheader, keepvcf, workdir, verbose)
            end
            msg = string("save parsed genofile: ", outfile)
            printconsole(logio,verbose,msg)
        else
            # isa(outstem, IO)
            parsevcf_io(inio, outstem, logio, pedinfo, isfounderinbred,
                formatpriority,isdelmultiallelic,commentstring, nheader, keepvcf, workdir, verbose)
        end
    end
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, logio, tused,"parvcf"; verbose,delim="-")
    return
end

function parsevcf_io(inio::IO,outio::IO,logio::Union{Nothing, IO},
    pedinfo::Union{Integer,AbstractString}, isfounderinbred::Bool,
    formatpriority::AbstractVector,
    isdelmultiallelic::Bool,    
    commentstring::AbstractString,
    nheader::Integer, 
    keepvcf::Bool,    
    workdir::AbstractString, verbose::Bool)
    filedate = string(commentstring*"fileDate=",replace(string(Date(now())),"-"=>""))
    write(outio, filedate,"\n")
    write(outio, commentstring*"source= modified by RABBIT/parsevcf\n")
    for _ in 1:nheader-1
        line = readline(inio,keep=true)        
        write(outio, line)
    end
    # parse title row
    titlerow = split(readline(inio,keep=false),"\t")
    magicped, fcols, offcols,newtitlerow = parse_titlerow(titlerow, pedinfo;
        keepvcf, commentstring,workdir,logio)
    write(outio,newtitlerow,"\n")
    # parse rest rows
    msg = string("begin, #founders=",length(fcols), ", #offspring=",length(offcols))
    printconsole(logio,verbose,msg)
    fgeno = Vector(undef, length(fcols))
    offgeno = Vector(undef, length(offcols))
    keys = Symbol.(formatpriority)
    values = 1:length(formatpriority)
    prioritytuple = (; zip(keys, values)...)
    missingset = [".","./.", ".|."]
    founder_nhetero = zeros(Int, length(fcols))
    nmarker = 0
    nmultia = 0
    outdelim = keepvcf ? '\t' : ','
    startt = time()
    while !eof(inio)        
        # rem(nmarker, 1000) == 0 && (startt = time())
        nmarker += 1
        rowgeno = split(readline(inio,keep=false),"\t")
        if length(rowgeno) <= 9
            @error string("#columns < 9 for rowgeno=",rowgeno, " at marker index=", nmarker)
        end
        # filter for multiallelic
        alt = rowgeno[5] # col5=alternative                
        ismultia =length(split(alt,",")) > 1                                
        if ismultia 
            nmultia += 1
            isdelmultiallelic && continue
        end          # 
        format = Symbol.(split(rowgeno[9],":"))
        formatcode = [get(prioritytuple,i,-1) for i in format] # -1 denotes missing fromat        
        founderformat = parse_rowgeno!(fgeno,rowgeno[fcols],
            formatcode, formatpriority, missingset,keepvcf)
        if isfounderinbred && (!keepvcf)
            ishetero, founderformat = to_GT_haplo!(fgeno,founderformat)
            if !ismissing(ishetero) 
                founder_nhetero .+= ishetero
            end
        end        
        offspringformat = parse_rowgeno!(offgeno, rowgeno[offcols],
            formatcode, formatpriority, missingset,keepvcf)                
        # if nmarker==1
        #     println("[rowgeno[offcols],formatcode, formatpriority, missingset,keepvcf]=",
        #         [rowgeno[offcols],formatcode, formatpriority, missingset,keepvcf])
        #     println("offgeno=",offgeno, ",offspringformat=",offspringformat)
        #     error("xx")
        # end
        if keepvcf
            if founderformat == offspringformat
                newformat = founderformat
            else
                newformat = string(founderformat,":", offspringformat)
            end
            # vcf column ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL","FILTER", "INFO", "FORMAT"]
            res = vcat(rowgeno[1:8], [newformat], fgeno, offgeno)
        else
            # keep information of some vcf columns            
            # parse vcf info column
            info = rowgeno[8]        
            newinfo, linkagegroup, poscm, foundererror, offspringerror, baseerror, allelicbias, allelicoverdispersion, allelicdropout = parse_vcf_info(info)
            idls = ["REF", "ALT", "QUAL","FILTER"]
            for i in eachindex(idls)
                if rowgeno[3+i] != "."
                    newinfo *= string(";",idls[i],"=",replace(rowgeno[3+i],","=>"|"))
                end
            end
            newinfo = strip(newinfo,';')    
            isempty(newinfo) && (newinfo = "NA")
            physchrom,physposbp,snpid = [i== "." ? "NA" : i for i in rowgeno[1:3]]            
            res = vcat([snpid,linkagegroup,poscm,physchrom,physposbp,newinfo], [founderformat, offspringformat],
                [foundererror, offspringerror,baseerror,allelicbias, allelicoverdispersion,allelicdropout], fgeno, offgeno)
            0
        end
        write(outio,join(res,outdelim),"\n")
        if rem(nmarker, 1000) == 0
            msg = string("marker index = ", nmarker, ", tused = ", 
                round(time()-startt,digits=1),"s")
            printconsole(logio,verbose,msg)
        end
    end
    nhetero = sum(founder_nhetero)
    if nhetero == 0
        if isfounderinbred
            msg = "no heterozygous founder genotypes detected from vcf to csv"
            printconsole(logio,verbose,msg)
        end
    else
        msg = string(nhetero, " out of ", length(fcols)*nmarker, " founder heterozygous genotypes are set to missing")
        @warn msg
        printconsole(logio,false,"Warning: "*msg)
        founderid = magicped.founderinfo[!,:individual]
        df = DataFrame(founder=founderid, nmarker=nmarker, nhetero=founder_nhetero)
        filter!(row->row.nhetero>0, df)
        sort!(df,:nhetero; rev=true)
        isa(logio, IO) && CSV.write(logio, df; header=true, append=true)
        msg = string(size(df,1), " out of ", length(founderid), " founders with  heterozygous genotypes")
        printconsole(logio,false,msg)
        println(msg, "\n", df)
    end
    msg = string("end, #founders=",length(fcols), ", #offspring=",length(offcols), ", #markers=",nmarker, 
        ",#multiallelic=",nmultia, isdelmultiallelic ? "(del)" : "(keep)")
    printconsole(logio,verbose,msg)
end

function parse_vcf_info(info::AbstractString)    
    newinfo = ""
    linkagegroup = "NA"          
    poscm = "NA"
    foundererror = "NA"
    offspringerror = "NA"
    baseerror = "NA"
    allelicbias = "NA"            
    allelicoverdispersion = "NA"            
    allelicdropout = "NA"
    if info != "."
        infols = split.(split(info,";"),"=")
        otherinfo = infols[length.(infols) .== 1]
        if !isempty(otherinfo)
            newinfo *=";"*join(reduce(vcat,otherinfo),";")
        end
        infols = infols[length.(infols) .== 2]
        infols = [strip.(i) for i in infols]                
        infodict = OrderedDict([Pair(i...) for i in infols])
        if haskey(infodict,"LINKAGEGROUP")
            linkagegroup = infodict["LINKAGEGROUP"]
            delete!(infodict,"LINKAGEGROUP")
        end
        if haskey(infodict,"POSCM")
            poscm = infodict["POSCM"]
            delete!(infodict,"POSCM")
        end
        if haskey(infodict,"FOUNDERERROR")
            foundererror = infodict["FOUNDERERROR"]
            delete!(infodict,"FOUNDERERROR")
        end
        if haskey(infodict,"OFFSPRINGERROR")
            offspringerror = infodict["OFFSPRINGERROR"]
            delete!(infodict,"OFFSPRINGERROR")
        end
        if haskey(infodict,"BASEERROR")
            baseerror = infodict["BASEERROR"]
            delete!(infodict,"BASEERROR")
        end
        if haskey(infodict,"ALLELICBIAS")
            allelicbias = infodict["ALLELICBIAS"]
            delete!(infodict,"ALLELICBIAS")
        end
        if haskey(infodict,"ALLELICOVERDISPERSION")
            allelicoverdispersion = infodict["ALLELICOVERDISPERSION"]
            delete!(infodict,"ALLELICOVERDISPERSION")
        end
        if haskey(infodict,"ALLELEDROPOUT")
            allelicdropout = infodict["ALLELEDROPOUT"]
            delete!(infodict,"ALLELEDROPOUT")
        end
        if !isempty(infodict)
            newinfo *=";"*join(map(i->string(i[1],"=",replace(i[2],","=>"|")),collect(infodict)),";")
        end
    end
    newinfo = strip(newinfo,';')    
    [newinfo, linkagegroup, poscm, foundererror, offspringerror, baseerror, allelicbias, allelicoverdispersion, allelicdropout]
end

function parse_titlerow(titlerow::AbstractVector,pedinfo::Union{Integer,AbstractString};
    keepvcf::Bool=true,
    commentstring::AbstractString='#',
    workdir::AbstractString=pwd(),    
    logio::Union{Nothing, IO})
    if isa(pedinfo, AbstractString)
        if last(splitext(pedinfo))==".csv"
            magicped = readmagicped(pedinfo; commentstring,workdir)
        else
            designinfo = parsedesign(pedinfo)
            nf = getnfounder(designinfo)
            length(titlerow) >= nf +9 || @error string("#columns(", length(titlerow),") < #founders(",nf, ") + 9")
            magicped = formmagicped!(designinfo,titlerow[10:end])
        end        
    else
        nf = pedinfo
        nf >=0 || @error string("Negative #founders=",nf)
        designinfo = parsedesign(string("nfounder=",nf))
        length(titlerow) >= nf +9 || @error string("#columns(", length(titlerow),") < #founders(",nf, ") + 9")        
        magicped = formmagicped!(designinfo,titlerow[10:end])
    end
    founderid = magicped.founderinfo[!,:individual]
    offspringid = magicped.offspringinfo[!,:individual]
    dict = Dict(titlerow .=> 1:length(titlerow))    
    fcols = [get(dict,i,nothing) for i in founderid]
    b = isnothing.(fcols)
    if any(b)
        msg = string(sum(b)," founders in pedfile but not in genofie: ",founderid[b])
        @warn msg
        printconsole(logio, false, "Warning: "*msg)
        fcols = fcols[.!b]
        founderid = founderid[.!b]
    end
    offcols = [if occursin(r"_virtualoffspring$",i)
        i2 = replace(i, r"_virtualoffspring$" => "")
        get(dict, i, get(dict,i2,nothing))
    else
        get(dict,i, nothing)
    end for i in offspringid]
    b = isnothing.(offcols)
    if any(b)        
        msg = string(sum(b)," offspring in pedfile but not in genofie: ",offspringid[b])
        @warn msg
        printconsole(logio, false, "Warning: "*msg)
        offcols = offcols[.!b]
        offspringid = offspringid[.!b]
    end
    if keepvcf
        res = vcat(titlerow[1:9], founderid, offspringid)
    else
        res = ["marker", "linkagegroup","poscm","physchrom", "physposbp", "info", "founderformat","offspringformat",
            "foundererror","offspringerror","baseerror","allelicbias","allelicoverdispersion","allelicdropout"]
        res = vcat(res,founderid, offspringid)
    end
    delim = keepvcf ? '\t' : ','    
    newtitlerow = join(res,delim)
    magicped, fcols, offcols,newtitlerow
end

function to_GT_haplo!(geno::AbstractVector,format::AbstractString)
    newformat = "GT_haplo"
    if format == "NA"
        ishetero = falses(length(geno))
        geno .= "N"
    elseif format == "GT_unphased"
        ishetero = falses(length(geno))
        for i in eachindex(geno,ishetero)
            gg = split(geno[i],"")
            if isempty(gg) 
                geno[i] = "N"
            elseif length(gg) == 1
                geno[i] = only(gg)
            elseif length(gg) == 2
                if allequal(gg)
                    geno[i] = first(gg)
                else
                    geno[i] = "N"
                    ishetero[i] = true
                end
            else
                @error string("unexpected genotype =",geno[i])
            end
        end
    elseif format == "GT_phased"
        ishetero = falses(length(geno))
        for i in eachindex(geno,ishetero)
            gg = split(geno[i],"|")
            if isempty(gg) 
                geno[i] = "N"
            elseif length(gg) == 1
                geno[i] = only(gg)
            elseif length(gg) == 2
                if allequal(gg)
                    geno[i] = first(gg)
                else
                    geno[i] = "N"
                    ishetero[i] = true
                end
            else
                @error string("unexpected genotype =",geno[i])
            end
        end
    elseif format in ["AD","GP"]
        # geno .= parse_AD2haplo.(geno)
        # ishetero = geno .== "-1"
        # geno[ishetero] .= "N"
        ishetero = missing
        newformat = format
    else
        @error string("TODO for format: ",format, ", geno=",geno)
    end
    ishetero, newformat
end

function parse_AD2haplo(g::AbstractString)
    # c1,c2 was transformed to c1&c2 by parse_vcfcell in  parse_rowgeno!
    if occursin(r"^0&[1-9][0-9]{0,}", g)
        "2"
    elseif occursin(r"^[1-9][0-9]{0,}&0$", g)
        "1"
    elseif occursin(r"^0&0$", g)
        "N"
    elseif occursin(r"^[1-9][0-9]{0,}&[1-9][0-9]{0,}$", g)
        "-1"
    elseif g == "NA"
        "N"
    else
        @error string("unknown AD=", g)
        "N"
    end
end


function parse_rowgeno!(resgeno::AbstractVector,
    rowgeno::AbstractVector,formatcode::AbstractVector,
    formatpriority::AbstractVector, missingset::AbstractVector,keepvcf::Bool)
    ncell = length(rowgeno)
    formatls = Vector(undef,ncell)
    if all(formatcode .< 0) #-1 denotes missing fromat
        resgeno .= missing
        formatls .= missing
    else
        ThreadsX.foreach(eachindex(resgeno,formatls,rowgeno)) do i
            @inbounds resgeno[i], formatls[i] = parse_vcfcell(rowgeno[i], formatcode,
                formatpriority, missingset,keepvcf)
        end
    end
    commonformat, missval = parse_format(formatls,formatpriority; isvcfformat = keepvcf)       
    resgeno[ismissing.(resgeno)] .= missval
    if ismissing(commonformat)
        if keepvcf
            commonformat = "GT"
            resgeno .= "."
        else
            commonformat = "GT_unphased"
            resgeno .= "NN"
        end
    else
        b =[i !== commonformat for i in formatls]
        resgeno[b] .= missval
    end
    commonformat
end

function parse_format(formatls::AbstractVector,formatpriority::AbstractVector; isvcfformat::Bool)
    ls = skipmissing(formatls)
    if isempty(ls)
        (missing,missing)
    else
        # commonformat = commonest(ls)               
        # assuming that only frist 2 characters of formats  are retained in formatpriority
        formatls2 = unique(ls)
        posls = [findfirst(formatpriority .== i[1:2]) for i in formatls2]
        commonformat =  formatls2[argmin(posls)]
        if isvcfformat
            missval = "."
        else
            if commonformat == "GT_unphased"
                missval = "NN"
            elseif commonformat == "GT_phased"
                missval = "N|N"
            else
                missval = "NA"
            end
        end
        commonformat,missval
    end
end

function parse_vcfcell(cellgeno::AbstractString,formatcode::AbstractVector,
    formatpriority::AbstractVector, missingset::AbstractVector,keepvcf::Bool)
    n = length(formatcode)
    geno = split(cellgeno,":")
    ng = length(geno)
    if ng < n
        geno = vcat(geno, repeat(["."],n-ng))
    elseif ng > n
        geno = geno[1:n]
    end
    b = [(geno[i] in missingset) || formatcode[i] < 0 for i in 1:n]  # -1 of formatcode denotes missing fromat
    if all(b)
        genoset = geno[formatcode .>= 1]
        if in("./.", genoset)
            g = "./."
            f = "GT" 
        elseif in(".|.", genoset)
            g = ".|."
            f = "GT" 
        else
            g = missing
            f = missing
        end
    else
        # the less the formatcode element, the higher the priority
        upbound = length(formatpriority)+10
        i = argmin(@. formatcode + upbound * b)             
        g = geno[i]
        f = formatpriority[formatcode[i]]
    end
    if keepvcf || ismissing(f)
        g, f
    else
        if f == "GT"
            if occursin("|", g)
                f2 = "GT_phased"
            elseif occursin("/", g)
                f2 = "GT_unphased"
            else
                error(string("unknown GT geno: ",g, ",cellgeno=",cellgeno))
            end
            # refenrence allele =1, alternative alelles = 2
            dict=Dict("0"=>"1","1"=>"2", "."=>"N","/"=>"","|"=>"|")
            g2 = join([get(dict,i,"2") for i=split(g,"")])
        else
            g2 = replace(g, "," => "&")
            f2 = f
        end
        g2, f2
    end
end


function tocsvgeno(genofile::AbstractString,pedinfo::Union{Integer,AbstractString};
    outstem::AbstractString = first(split_allext(genofile)),
    outext::AbstractString = ".csv.gz",
    isfounderinbred::Bool=true,
    formatpriority::AbstractVector=["AD","GT"],
    commentstring::AbstractString="##",
    workdir::AbstractString=pwd())
    last(split_allext(genofile)) in [".vcf",".vcf.gz"] || @error string("inpu genofile ext must .vcf or .vcf.gz")    
    parsevcf(genofile, pedinfo;
            outstem, outext,isfounderinbred,
            formatpriority, commentstring, workdir, logfile = nothing, verbose=false)
    outstem*outext
end

function tovcfgeno(genofile::AbstractString,pedinfo::Union{Integer,AbstractString};
    outstem::AbstractString = first(split_allext(genofile)),
    outext::AbstractString = ".vcf.gz",
    isfounderinbred::Bool=true,
    commentstring::AbstractString="##",
    workdir::AbstractString=pwd())
    last(split_allext(genofile)) in [".csv",".csv.gz"] || @error string("input genofile ext must .csv or .csv.gz")
    magicgeno = formmagicgeno(genofile, pedinfo;
        isfounderinbred,isphysmap=false,
        commentstring, missingstring=["NA","missing"], workdir)
    outfile = outstem*outext    
    savegenodata(outfile,magicgeno;
        delim='\t',commentstring,workdir)
    outfile
end