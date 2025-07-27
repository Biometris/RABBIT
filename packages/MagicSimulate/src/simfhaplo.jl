"""
    simfhaplo(; kwargs...)

simulate founder haplotypes. 

# Keyword arguments

`nsnp`: number of markers.

`nparent`: number of parents.

`missingstring="NA"`: string representing missing value. 

`chrlen::AbstractVector = 100 * ones(5)`:  specify genetic length of each chromosome. 

`isfounderinbred::Bool=true`: if true, founders are inbred, and otherwise outbred.

`outfile = "sim_fhaplo.vcf.gz"`: output filename. 

`workdir::AbstractString = pwd()` specifies the working directory.

"""
function simfhaplo(;
    nsnp::Integer,
    nparent::Integer,
    chrlen::AbstractVector = 100 * ones(5),
    isfounderinbred= true,
    missingstring="NA",
    outfile = "sim_fhaplo.vcf.gz",
    workdir::AbstractString = pwd())
    fileext = last(MagicBase.split_allext(outfile))
    fileext in [".csv",".csv.gz",".vcf",".vcf.gz"] || @error string("outfile ext must be in [.csv,.csv.gz,.vcf,.vcf.gz]")
    isvcf = fileext in [".vcf",".vcf.gz"]
    myopen = fileext == ".csv.gz" ? GZip.open : open
    if isvcf
        # temp csv file will be transformed into output vcf file
        outfile2 = tempname(workdir;cleanup=true)*".csv"
    else
        outfile2 = MagicBase.getabsfile(workdir,outfile)
    end
    try
        myopen(outfile2, "w") do io
            commentstring = "##"        
            write(io, commentstring*"fileformat=",fileext,"\n")
            filedate = string(commentstring*"fileDate=",replace(string(Date(now())),"-"=>""))
            write(io, filedate,"\n")
            write(io, commentstring*"source=MagicSimulate\n")
            nsnpls = round.(Int,nsnp .* (chrlen ./ sum(chrlen)))
            nsnpls[1] += nsnp - sum(nsnpls)
            for chr in 1:length(chrlen)
                nsnp = nsnpls[chr]
                d = rand(Exponential(1),nsnp-1)
                pos = accumulate(+,d)
                pushfirst!(pos,0)
                snppos = round.(pos ./ (last(pos)/chrlen[chr]), digits=4)
                snp0 = sum(nsnpls[1:chr-1])
                snpid = [string("snp",i+snp0) for i in 1:nsnp]
                parentid = [string("P",i) for i in 1:nparent]
                founderformat = isfounderinbred ? "GT_haplo" : "GT_phased"
                markermap = DataFrame("marker" => snpid,
                    "linkagegroup" => chr, 
                    "poscm" => snppos,
                    "physchrom" => missing,
                    "physposbp" => missing,
                    "info" => missing,
                    "founderformat" => founderformat,
                    "offspringformat" => missing,
                    "foundererror"=>missing,
                    "offspringerror"=>missing,
                    "baseerror"=>missing,
                    "allelicbias"=>missing,
                    "allelicoverdispersion"=>missing,
                    "allelicdropout"=>missing,
                )
                nfgl = isfounderinbred ? nparent : 2*nparent
                haplo = ones(Int, nfgl,nsnp)
                for j in 1:nsnp
                    nallele2= sample(1:nfgl-1)
                    ii = sample(1:nfgl, nallele2, replace = false)
                    haplo[ii,j] .= 2
                end
                haplo = haplo'
                if isfounderinbred
                    haplodf = DataFrame(string.(haplo),parentid)
                else
                    haplo = haplo[:,1:2:end] .* 10 + haplo[:,2:2:end]
                    haplodf = DataFrame(join.(digits.(haplo),"|"),parentid)
                end
                df = hcat(markermap,haplodf)                
                CSV.write(io, df; delim=',', missingstring, header=(chr==1),append=true)
            end
        end
        if isvcf
            resultfile = MagicBase.getabsfile(workdir,outfile)
            outstem,outext = MagicBase.split_allext(resultfile)
            MagicBase.tovcfgeno(outfile2,nparent; outstem, outext,isfounderinbred)
        else
            resultfile = outfile2
        end
        resultfile
    finally
        isvcf && rm(outfile2)
    end
end
