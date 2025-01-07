"""
    hapmap2vcf(hapmapfile; keyargs...)

convert from a hapmap genofile into a vcf genofile. 

# Positional arguments

`hapmapfile::AbstractString`: hapmap genofile. 

# Keyword arguments

`delmultiallelic::Bool = true,`: if true, delete markers with #alleles >= 3. 

`delim= ","`: text delimiter. 

`missingallele::AbstractString="-"`: string for missing allele

`outstem::AbstractString="outstem"`: stem of output filename.

`workdir::AbstractString=pwd()`: directory for reading and saving files.

`logfile::Union{Nothing,AbstractString,IO}= outstem*"_arrayfile2vcf.log"`: log file or IO for writing log. If it is nothing, no log file.

`verbose::Bool=true`: if true, print details on the stdout.

"""
function hapmap2vcf(hapmapfile::AbstractString;     
    delmultiallelic::Bool = true,
    missingallele::AbstractString="-",
    delim = '\t',    
    outstem::AbstractString  = "outstem", 
    workdir::AbstractString = pwd(),
    logfile::Union{Nothing,AbstractString,IO}= outstem*"_hapmap2vcf.log",
    verbose::Bool=true)
    starttime = time()
    logio = MagicBase.set_logfile_begin(logfile, workdir, "hapmap2vcf"; verbose)       
    msg = string("list of args: \n",
        "hapmapfile = ", hapmapfile, "\n",                
        "delim = ", string(delim), "\n",        
        "workdir = ",workdir,"\n",        
        "outstem = ", isnothing(outstem) ? "no output files" : outstem,"\n",        
        "logfile = ",logfile,"\n",
        "verbose = ",verbose)
    printconsole(logio,verbose,msg)        
    hapmapfile2 = getabsfile(workdir,hapmapfile)    
    occursin(r".hmp.txt$",hapmapfile2) || @error string("file extension  is not .hmp.txt for hapmapfile=",hapmapfile)
    nmarker = 0
    nmulti = 0
    nsample = [0]
    outfile = getabsfile(workdir,outstem*"_geno.vcf.gz")        
    open(hapmapfile2,"r") do inio    
        GZip.open(outfile, "w") do outio    
            title = split(readline(inio;keep=false),delim)
            length(title) >= 11 || @error string("#columns in hapmapfile = ",length(title), " < 11, title=",title)
            hapmapcols = ["rs#","alleles","chrom","pos","strand","assembly#","center",    
                "protLSID","assayLSID","panelLSID","QCcode"] 
            if title[1:11] != hapmapcols
                @error string("the first 11 columns: ", title[1:11], 
                    " not equal to ", join(hapmapcols,"",))
            end
            samples = title[12:end]
            nsample[1] = length(samples) 
            msg = string(nsample[1], " genotyped individuals (print upto 10): ", join(samples[1:min(10,length(samples))],","))
            printconsole(logio, verbose, msg)
            title =vcat(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL","FILTER", "INFO", "FORMAT"],samples)
            write(outio, join(title,'\t'),"\n")     
            nwarnmax = 10
            nwarn = 0            
            while !eof(inio)    
                nmarker += 1
                rowgeno = split(readline(inio,keep=false),delim)
                markerid,chrom, physposbp = rowgeno[[1,3,4]]
                # calculate alleles with decreasing frequencies                
                genostr = replace(string(rowgeno[12:end]...),missingallele=>"")
                alleles = string.(unique(genostr))
                # https://www.bioinformatics.org/sms/iupac.html                
                issubset(alleles,["A","T","G","C","R","Y","S","W","K","M","B","D","H","V","N","."]) || @error string("unknown genotype codes =",alleles)
                intersect!(alleles,["A","T","G","C"])
                freqls = [count(i,genostr) for i in alleles]
                oo = sortperm(freqls,rev=true)
                alleles .= alleles[oo]                
                if isempty(alleles)
                    # VCFv4.3 requires that REF must be one of A, T, G, C, N
                    ref = "N"
                    alt = "."                    
                elseif length(alleles) == 1                    
                    ref = only(alleles)
                    alt = "."
                else
                    if length(alleles) > 2
                        nmulti += 1
                        delmultiallelic && continue
                    end                    
                    ref = first(alleles)
                    alt = join(alleles[2:end],"")
                end
                info = ""
                for j in 5:11
                    if rowgeno[j] != "NA" 
                        if isempty(info) 
                            info *= string(hapmapcols[j],"=",rowgeno[j])
                        else
                            info *= string(";", hapmapcols[j],"=",rowgeno[j])
                        end
                    end
                end
                isempty(info) && (info = ".")                
                col19 = [chrom,physposbp,markerid,ref,alt,".",".",info,"GT"]                      
                allele_rule = alleles .=> string.(0:(length(alleles)-1)) 
                for j in 12:length(rowgeno)
                    if occursin(r"[R|Y|S|W|K|M|B|D|H|V]",rowgeno[j])                        
                        if nwarn <= nwarnmax
                            @warn string("code R|Y|S|W|K|M|B|D|H|V is set missing")      
                            nwarn += 1     
                        elseif nwarn+1 == nwarnmax
                            @warn string("stop warning: code R|Y|S|W|K|M|B|D|H|V is set missing")      
                        else
                            0
                        end
                    end                    
                    g = split(replace(rowgeno[j],r"[R|Y|S|W|K|M|B|D|H|V|N|.]"=>".",allele_rule...),"")      
                    rowgeno[j] = join(g,"/")
                end
                write(outio, join(col19,"\t"), "\t", join(rowgeno[12:end],"\t"),"\n")
            end
        end
    end
    msg = string("#individuals=",nsample[1],", #markers=",nmarker)
    if nmulti > 0
        msg *= string(", #markers_multiallelic=",nmulti, delmultiallelic ? "(deleted)" : "" )
    end
    printconsole(logio,verbose,msg)    
    msg = string("ouput vcffile: ", outfile)
    printconsole(logio,verbose,msg)    
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, logio, tused,"hapmap2vcf"; verbose)        
    outfile
end
