"""
    magicsimulate(fhaplofile, pedinfo; kwargs...)

simulate genotypic data from founder haplotypes in `fhaplofile` and pedigree
information in `pedinfo`.

# Positional arguments

`fhaplofile::AbstractString` specifies file for founder haplotypes including marker map.
The fhaplofile extension must be ".vcf" or ".csv".

`pedinfo::Union{AbstractString,MateScheme}` specifies pedigree information via a matescheme
(e.g. MateScheme(8,["Pairing","Selfing"],[3,6])), 
or via a pedigree file (if pedinfo[end-3:end]==".csv"), 
or via a string designcode (e.g. "8ril-self6").

# Keyword arguments

`isfounderinbred::Bool=true` if true, the founders are inbred. For inbred founders,
the heterozygous genotypes are set to missing. For outbred founders, genotypes must
be phased.

`popsize::Union{Nothing,Integer}=200` specifies the population size, valid only if
the `pedinfo` is specified via a designcode or matescheme.

`error_randallele::Union{Nothing,Real}=nothing` specifices genotyping error model. 
The error model follows uniform allele model with probability error_randallele,  and it follows
uniform genotype model with probability 1-error_randallele. If it is nothing, error_randallele 
is given by the interally estimated non-ibd probability. 

`foundererror::Distribution=Beta(1,199)` specifies that the probability distribution
of genotyping error rate at a marker in founders. At a given marker, founders have the same error rate. 

`offspringerror::Distribution=Beta(1,199)` specifies that the probability distribution
of genotyping error rate at a marker in offspring. At a given marker, offspring have the same error rate. 

`foundermiss::Distribution=Beta(1,9)` specifies that the probability distribution
of the fraction of missing genotypes at a marker in founders.

`offspringmiss::Distribution=Beta(1,9)` specifies that the probability distribution
of the fraction of missing genotypes at a marker in offspring.

`seqfrac::Real=0.0` specifies the fraction of markers being genotyped by sequencing;
the rest markers are genotyped by SNP array.

`baseerror::Distribution = Beta(1,199)` specifies the probability diestribution of
sequencing error rate among markers. 

`allelicbias::Distribution = Beta(10,10)` specifies the probability distribution of
the mean sequencing allelic balance among markers. 

`allelicoverdispersion::Distribution = Exponential(0.05)` specifies overdispersion parameter for the probability distribution of
allelicbias among individuals at a marker. 

`seqdepth::Distribution = Gamma(2, 5)` specifies the probability distribution
of the mean read depth at a marker. 

`seqdepth_overdispersion::Distribution = Gamma(1,1)` specifies the probability distribution
of over-dispersion of read depths among individuals at a marker. 
Given the mean depth lam at a marker, the read depth of an individual
follows a NegativeBionomial(r,p), such that mean lam = r(1-p)/p, and variance = r(1-p)/p^2 = r(1+lam/r) = r(1+seqdepth_overdispersion)
where seqdepth_overdispersion = lam/r. seqdepth_overdispersion = 0 denotes no over-dispersion. 

`isobligate::Bool=false` specifies whether there must be at least one crossover event

`interference::Integer=0` specifies the strength of chiasma interference. By default,
no chiasma interference. The recombination breakpoins are obtained by taking every
(1+interference) points that follow a Poisson distribution along chromosome.

`ispheno::Bool=false` specifies whether to simulate phenotypes.

`pheno_nqtl::Integer=1` specifies the number of QTLs in simulating phenotypes.

`pheno_h2::Real= 0.5` specifies heritablity in simulating phenotypes.

`select_nqtl::Integer=1` specifies the number of QTLs in simulating trait for artifical selection.

`select_prop::Real = 1.0` specifies the proportion of zygotes selected in artifical selection.
By default, no artifical selection.

`outstem::Union{Nothing,AbstractString}="outstem"` specifies the stem of output filenames.

`workdir::AbstractString = pwd()` specifies the working directory.

`verbose::Bool=true`: if true, print details on the stdout.

# Examples
```julia-repl
julia> magicsimulate("fhaplo.vcf.gz","ped.csv")
julia> magicsimulate("fhaplo.vcf.gz","8ril-self6"; popsize=800)
```
"""
function magicsimulate(fhaplofile::AbstractString,
    pedinfo::Union{AbstractString,MateScheme};
    isfounderinbred::Bool=true,
    popsize::Union{Nothing,Integer}=200,
    foundermiss::Distribution=Beta(1,9),
    offspringmiss::Distribution=Beta(1,9),
    error_randallele::Union{Nothing,Real}=nothing,
    foundererror::Distribution=Beta(1,199),
    offspringerror::Distribution=Beta(1,199),    
    seqfrac::Real=0.0,    
    baseerror::Distribution = Beta(1,199),
    allelicbias::Distribution = Beta(10,10),
    allelicoverdispersion::Distribution = Exponential(0.05),
    allelicdropout::Distribution = Uniform(0,1e-10),
    seqdepth::Distribution = Gamma(2, 5),
    seqdepth_overdispersion::Distribution = Gamma(1,1),
    isobligate::Bool=false,
    interference::Integer=0,
    ispheno::Bool=false,
    pheno_nqtl::Integer=1,    
    pheno_h2::Real= 0.5,
    select_nqtl::Integer=50,
    select_dom::Real=0, 
    select_prop::Real = 1.0, # 1.0 === no artifical selection
    nplot_subpop::Integer=10,
    commentstring::AbstractString="##",
    missingstring::AbstractString="NA",
    outstem::Union{Nothing,AbstractString}="outstem",
    outext::AbstractString=".vcf.gz",
    workdir::AbstractString = pwd(),
    logfile::Union{Nothing,AbstractString,IO}= (isnothing(outstem) ? nothing : string(outstem,"_magicsimulate.log")),
    verbose::Bool=false)
    starttime = time()
    io = MagicBase.set_logfile_begin(logfile, workdir, "magicsimulate"; verbose)    
    if isa(logfile, AbstractString)
        msg=string("Julia version: ",VERSION)
        MagicBase.printconsole(io,verbose,msg)
        MagicBase.printpkgst(io,verbose,"MagicSimulate")
    end
    MagicBase.check_infiles(fhaplofile,pedinfo; isbreedped=false, ismagicsimulate=true, io, commentstring,workdir,verbose)
    check_sim_arg(error_randallele, foundererror, offspringerror, foundermiss, offspringmiss,
        seqfrac, baseerror, allelicbias, allelicoverdispersion, allelicdropout, seqdepth, seqdepth_overdispersion, 
        popsize, interference, pheno_nqtl, pheno_h2,
        select_nqtl, select_dom, select_prop,outext)
    msg = string("list of args: \n",
        "fhaplofile = ", fhaplofile, "\n",
        "pedinfo = ", pedinfo, "\n",
        "isfounderinbred = ", isfounderinbred,"\n",
        "popsize = ", popsize, "\n",
        "error_randallele = ", error_randallele, "\n",
        "foundererror = ", foundererror, "\n",
        "offspringerror = ", offspringerror, "\n",
        "foundermiss = ", foundermiss, "\n",
        "offspringmiss = ", offspringmiss, "\n",            
        "seqfrac = ", seqfrac,"\n",
        "baseerror = ", baseerror, "\n",
        "allelicbias = ", allelicbias, "\n",
        "allelicoverdispersion = ", allelicoverdispersion, "\n",
        "allelicdropout = ", allelicdropout, "\n",
        "seqdepth = ", seqdepth,"\n",
        "seqdepth_overdispersion = ", seqdepth_overdispersion,"\n",        
        "isobligate = ", isobligate, "\n",
        "interference = ", interference, "\n",
        "ispheno = ", ispheno, "\n",
        "pheno_nqtl = ", pheno_nqtl,"\n",
        "pheno_h2 = ", pheno_h2,"\n",
        "select_nqtl = ", select_nqtl,"\n",
        "select_dom = ", select_dom,"\n",
        "select_prop = ", select_prop,"\n",
        "nplot_subpop = ", nplot_subpop,"\n",
        "commentstring = ", commentstring,"\n",
        "missingstring = ", missingstring,"\n",
        "workdir = ",workdir,"\n",
        "outstem = ", isnothing(outstem) ? "no output files" : outstem,"\n",        
        "outext = ",outext,"\n",
        "logfile = ",logfile,"\n",        
        "verbose = ",verbose)
    MagicBase.printconsole(io,verbose,msg)
    # 
    tused = @elapsed begin 
        tused = @elapsed magicped = pedinfo2magicped(pedinfo; popsize,isfounderinbred, workdir) 
        if isa(magicped.designinfo,Pedigree)
            magicped.offspringinfo[!,:isfglexch] .= false
        end
        nfounder = size(magicped.founderinfo,1)
        founderhaplo = MagicBase.formfhaplo(fhaplofile; 
            formatpriority = ["GT"], isphysmap=false, 
            missingstring = unique(vcat(missingstring, "missing")),
            commentstring, workdir,isfounderinbred)        
        nfounder_haplo = size(founderhaplo.magicped.founderinfo,1) 
        nfounder_haplo < nfounder  && @error string("#founders=", nfounder_haplo," in fhaplofile=", " < #founders=",nfounder, " in pedfile")
        if nfounder_haplo > nfounder 
            founderhaplo.magicped.founderinfo  = founderhaplo.magicped.founderinfo[1:nfounder,:]
            founderhaplo.foundergeno  .= [i[:,1:nfounder] for i in founderhaplo.foundergeno]
        end
        if !ispedfile(pedinfo)
            MagicBase.setfounderid!(magicped,founderhaplo.magicped.founderinfo[!,:individual])
        end
        MagicBase.rawgenoprob!(founderhaplo; targets = ["founders"], baseerror=0.005, isfounderinbred)
        MagicBase.rawgenocall!(founderhaplo; targets = ["founders"], callthreshold = 0.95, isfounderinbred)
        founderformat = unique(vcat([unique(i[!,:founderformat]) for i=founderhaplo.markermap]...))
        formatdiff = setdiff(founderformat, ["GT_haplo","GT_phased"])
        if !isempty(formatdiff)
            @error string("unexpected founderformat: ", formatdiff)
        end
    end    
    msg = string("tused=",round(tused,digits=2),"s in reading input files")
    MagicBase.printconsole(io,verbose,msg)
    #     
    tused = @elapsed contfgl, magicfgl = simfgl(founderhaplo,magicped;
        isfounderinbred,isobligate,interference,
        select_nqtl,select_prop)
    msg = string("tused=",round(tused,digits=2),"s in simulating inheritance pattern from pedinfo")
    MagicBase.printconsole(io,verbose,msg)
    # 
    tused = @elapsed truegeno = genedrop(magicfgl; isfounderinbred,isphased=true)
    msg = string("tused=",round(tused,digits=2),"s in simulating truegeno from inheritance pattern and fhaplo")
    MagicBase.printconsole(io,verbose,msg)
    # @info truegeno.markermap[1]
    # 
    tused = @elapsed obsgeno = simgeno(truegeno,magicfgl;error_randallele,foundererror,offspringerror,
        foundermiss,offspringmiss,seqfrac,baseerror,allelicbias,allelicoverdispersion,allelicdropout,
        seqdepth,seqdepth_overdispersion)
    msg = string("tused=",round(tused,digits=2),"s in simulating error and missing")
    MagicBase.printconsole(io,verbose,msg)
    # 
    if ispheno
        tused = @elapsed phenores = simpheno(contfgl,magicfgl; pheno_nqtl,pheno_h2)
        msg = string("tused=",round(tused,digits=2),"s in simulating phenotypes")
        MagicBase.printconsole(io,verbose,msg)
    else
        phenores = nothing
    end
    if !isnothing(outstem)
        startt =time()
        outstem *= "_magicsimulate"
        idls = ["geno","ped","truegeno","truefgl","truecontfgl"]
        outfiles = [string(outstem,"_",id,".csv.gz") for id=idls]
        delim = outext in [".vcf",".vcf.gz"] ? '\t' : ','
        outfiles[1] = string(outstem,"_geno",outext)
        outfiles[2] = string(outstem,"_ped.csv")
        savegenodata(outfiles[1],obsgeno; missingstring, workdir, delim)
        savemagicped(outfiles[2],magicped; workdir)
        savegenodata(outfiles[3],truegeno; missingstring, workdir)
        savegenodata(outfiles[4],magicfgl; missingstring, workdir)
        savegenodata(outfiles[5],contfgl; missingstring, workdir)        
        pedfile = string(outstem,"_ped.png")
        pedfile2 = MagicBase.getabsfile(workdir,pedfile)
        nondiv = allunique(magicped.offspringinfo[!,:member])
        isplotped = (!nondiv) || (length(magicped.designinfo.member)<1000)
        if isplotped
            try 
                plotmagicped(magicped;isfounderinbred,outfile=pedfile2)
                push!(outfiles,pedfile)
            catch
                @warn string("Could not plot magicped")
            end             
        end
        subpop2off = MagicBase.get_subpop2offspring(magicped; isindex=true)
        offls = sort(reduce(vcat,[length(val)<=nplot_subpop ? val : val[1:nplot_subpop] for (key, val) in subpop2off]))                           
        try 
            mosaicfiles = savemosaic(contfgl; 
                workdir,offspring=offls, 
                outstem=outstem*"_truecontfgl",io,verbose)             
            push!(outfiles,mosaicfiles...)
        catch err
            @warn string("Could not plot mosaic")
            @error err
        end    
        if ispheno
            fgl_qtl,geno_qtl,gadd_qtl,map_qtl,pheno_df = phenores
            # save pheno
            outfile = outstem*"_pheno.csv"
            push!(outfiles,outfile)
            CSV.write(MagicBase.getabsfile(workdir,outfile), pheno_df; missingstring, header=true)
            # save pheno truevalues
            outfile = outstem*"_truepheno.csv"
            push!(outfiles,outfile)
            open(MagicBase.getabsfile(workdir,outfile),"w") do outio
                write(outio,"MagicSimulate,map_qtl\n")
                CSV.write(outio, map_qtl; header=true,missingstring, append=true)
                resls = [fgl_qtl,geno_qtl,gadd_qtl]
                residls = ["fgl_qtl","geno_qtl","gadd_qtl"]
                for i in eachindex(resls)
                    write(outio,string("MagicSimulate,",residls[i],"\n"))
                    CSV.write(outio, resls[i]; header=true,missingstring, append=true)
                end
            end
        end
        msg = string("tused=",round(time()-startt,digits=2),"s in exporting")
        MagicBase.printconsole(io,verbose,msg)
        tused = round(time()-starttime,digits=2)
        MagicBase.set_logfile_end(logfile, io, tused,"magicsimulate"; verbose,delim="-")
        outfiles
    else
        tused = round(time()-starttime,digits=2)
        MagicBase.set_logfile_end(logfile, io, tused,"magicsimulate"; verbose,delim="-")
        (obsgeno=obsgeno, truegeno=truegeno, magicfgl=magicfgl, contfgl=contfgl,magicped=magicped, phenores=phenores)
    end
end

function check_sim_arg(
    error_randallele::Union{Nothing,Real},
    foundererror::Distribution,
    offspringerror::Distribution,
    foundermiss::Distribution,
    offspringmiss::Distribution,    
    seqfrac::Real,
    baseerror::Distribution,
    allelicbias::Distribution,
    allelicoverdispersion::Distribution, 
    allelicdropout::Distribution,
    seqdepth::Distribution,
    seqdepth_overdispersion::Distribution,
    popsize::Integer,    
    interference::Integer,
    pheno_nqtl::Integer,
    pheno_h2::Real,
    select_nqtl::Integer,
    select_dom::Real,
    select_prop::Real,
    outext::AbstractString)
    if !(isnothing(error_randallele) || 0<=error_randallele<=1)
        @error string("error_randallele=",error_randallele, "; it must be nothing or in [0,1]")
    end
    if !(minimum(foundererror) >=0 &&  maximum(foundererror) <=1)
        @error string("foundererror=",foundererror, "; its support is not a subset of [0,1]")
    end
    if !(minimum(offspringerror) >=0 &&  maximum(offspringerror) <=1)
        @error string("offspringerror=",offspringerror, "; its support is not a subset of [0,1]")
    end
    if !(minimum(foundermiss) >=0 &&  maximum(foundermiss) <=1)
        @error string("foundermiss=",foundermiss, "; its support is not a subset of [0,1]")
    end
    if !(minimum(offspringmiss) >=0 &&  maximum(offspringmiss) <=1)
        @error string("offspringmiss=",offspringmiss, "; its support is not a subset of [0,1]")
    end
    if !(0<=seqfrac<=1)
        @error string("seqfrac=", seqfrac,"; fraction of markers being sequenced is not in [0,1]")
    end
    if !(popsize>0)
        @error string("popsize=", popsize," is not positive")
    end
    if !(minimum(baseerror) >=0 &&  maximum(baseerror) <=1)
        @error string("baseerror=",baseerror, "; its support is not a subset of [0,1]")
    end
    if !(minimum(allelicbias) >=0 &&  maximum(allelicbias) <=1)
        @error string("allelicbias=",allelicbias, "; its support is not a subset of [0,1]")
    end
    if !(minimum(allelicoverdispersion) >=0)
        @error string("allelicoverdispersion=",allelicoverdispersion, "; its support is not a subset of [0,∞]")
    end
    if !(minimum(allelicdropout) >=0 &&  maximum(allelicdropout) <=1)
        @error string("allelicdropout=",allelicdropout, "; its support is not a subset of [0,1]")
    end
    if !(minimum(seqdepth) >= 0)
        @error string("seqdepth=", seqdepth,"; its support is not a subset of [0,∞)")
    end
    if !(minimum(seqdepth_overdispersion) >= 0)
        @error string("seqdepth_overdispersion=", seqdepth_overdispersion,"; its support is not a subset of [0,∞)")
    end    
    if !(interference>=0)
        @error string("interference=", interference," is not non-negative")
    end
    if !(pheno_nqtl>0)
        @error string("pheno_nqtl=", pheno_nqtl," is not positive")
    end
    if !(0<=pheno_h2<1)
        @error string("pheno_h2=", pheno_h2," is not [0,1)")
    end
    if !(select_nqtl>0)
        @error string("select_nqtl=", select_nqtl," is not positive")
    end
    if !(select_dom>=0)
        @error string("select_dom=", select_dom," is not non-negative")
    end
    if !(0<=select_prop<=1)
        @error string("select_prop=", select_prop," is not in [0,1]")
    end
    if !in(outext, [".vcf.gz",".csv.gz", ".vcf",".csv"])
        @error string("outext=",outext, ", not in [.vcf.gz, .csv.gz, .vcf, .csv]")
    end
    nothing
end

"""
    magicsimulate(pedinfo; kwargs...)

simulates ancestral blocks from pedinfo. 

# Positional arguments

`pedinfo::Union{AbstractString,MateScheme}` specifies pedigree information via a matescheme
(e.g. MateScheme(8,["Pairing","Selfing"],[3,6])), 
or via a pedigree file (if pedinfo[end-3:end]==".csv"), 
or via a string designcode (e.g. "8ril-self6").

# Keyword arguments

`isfounderinbred::Bool=true` if true, the founders are inbred. For inbred founders,
the heterozygous genotypes are set to missing. For outbred founders, genotypes must
be phased.

`popsize::Union{Nothing,Integer}=200` specifies the population size, valid only if
the `pedinfo` is specified via a designcode or matescheme.

`chrlen::AbstractVector= 100*ones(5)` specifies lengths (cM) for each chromosome.

`isobligate::Bool=false` specifies whether there must be at least one crossover event

`interference::Integer=0` specifies the strength of chiasma interference. By default,
no chiasma interference. The recombination breakpoins are obtained by taking every
(1+interference) points that follow a Poisson distribution along chromosome.

`outstem::Union{Nothing,AbstractString}="outstem"` specifies the stem of output filenames.

`workdir::AbstractString = pwd()` specifies the working directory.

`verbose::Bool=true`: if true, print details on the stdout.

# Outputs

| Output file                  | Description |
|:-----------------------------|:----------- |
|`outstem_ped.csv`        |  simulated pedigree file |
|`outstem_truecontfgl.csv`|  truevalues of continuous origin-genotypes |

Here fgl denotes founder genome labels, and origin-genotypes denote genotypes
with each fgl  being regarded as a distinct allele.

# Examples
```julia-repl
julia> magicsimulate("8ril-self6"; popsize=800)
```
"""
function magicsimulate(pedinfo::Union{AbstractString,MateScheme};
    isfounderinbred::Bool=true,
    popsize::Union{Nothing,Integer}=200,    
    chrlen::AbstractVector= 100*ones(5), # centiMorgan
    isobligate::Bool=false,
    interference::Integer=0,
    nplot_subpop::Integer=10,
    missingstring::AbstractString="NA",
    outstem::Union{Nothing,AbstractString}="outstem",
    workdir::AbstractString = pwd(),
    logfile::Union{Nothing,AbstractString,IO}= (isnothing(outstem) ? nothing : string(outstem,"_magicsimulate.log")),
    verbose::Bool=false)
    starttime = time()
    io = MagicBase.set_logfile_begin(logfile, workdir, "magicsimulate"; verbose)
    if isa(logfile, AbstractString)
        msg=string("Julia version: ",VERSION)
        MagicBase.printconsole(io,verbose,msg)
        MagicBase.printpkgst(io,verbose,"MagicSimulate")
    end
    msg = string("list of args: \n",       
        "pedinfo = ", pedinfo, "\n",
        "isfounderinbred = ", isfounderinbred,"\n",
        "popsize = ", popsize, "\n",        
        "chrlen = ", chrlen, "\n",
        "isobligate = ", isobligate, "\n",
        "interference = ", interference, "\n",
        "nplot_subpop = ", nplot_subpop,"\n",
        "missingstring = ", missingstring,"\n",
        "workdir = ",workdir,"\n",
        "outstem = ", isnothing(outstem) ? "no output files" : outstem,"\n",                
        "logfile = ",logfile,"\n",        
        "verbose = ",verbose)
    MagicBase.printconsole(io,verbose,msg)
    # 
    tused = @elapsed magicped = pedinfo2magicped(pedinfo; popsize,isfounderinbred,workdir)         
    msg = string("tused=",round(tused,digits=2),"s in reading pedigree info")
    MagicBase.printconsole(io,verbose,msg)
    # 
    tused = @elapsed begin 
        chrid = [string("chr",i) for i in eachindex(chrlen)]
        contfgl = simcontfgl(magicped; chrid,chrlen,
            isfounderinbred,isobligate,interference)        
    end
    msg = string("tused=",round(tused,digits=2),"s in simulating contfgl from pedinfo")
    MagicBase.printconsole(io,verbose,msg)
    # 
    if !isnothing(outstem)
        startt =time()
        outfiles = [outstem*"_ped.csv",outstem*"_truecontfgl.csv.gz"]        
        savemagicped(outfiles[1],magicped; workdir)
        savegenodata(outfiles[2],contfgl; missingstring, workdir)
        pedfile = string(outstem,"_ped.png")
        pedfile2 = MagicBase.getabsfile(workdir,pedfile)
        nondiv = allunique(magicped.offspringinfo[!,:member])
        isplotped = nplot_subpop > 0 && (!nondiv) || (length(magicped.designinfo.member)<1000)
        if isplotped
            try 
                plotmagicped(magicped; isfounderinbred, outfile=pedfile2)
                push!(outfiles,pedfile)
            catch
                @warn string("Could not plot magicped")
            end            
        end
        subpop2off = MagicBase.get_subpop2offspring(magicped; isindex=true)
        offls = sort(reduce(vcat,[length(val)<=nplot_subpop ? val : val[1:nplot_subpop] for (key, val) in subpop2off]))                           
        try 
            mosaicfiles = savemosaic(contfgl; 
                workdir,offspring=offls, 
                outstem=outstem*"_truecontfgl",io,verbose)             
            push!(outfiles,mosaicfiles...)
        catch
            @warn string("Could not plot mosaic")
        end
        msg = string("tused=",round(time()-startt,digits=2),"s in exporting")
        MagicBase.printconsole(io,verbose,msg)
        tused = round(time()-starttime,digits=2)
        MagicBase.set_logfile_end(logfile, io, tused,"magicsimulate"; verbose,delim="-")
        outfiles
    else
        tused = round(time()-starttime,digits=2)
        MagicBase.set_logfile_end(logfile, io, tused,"magicsimulate"; verbose,delim="-")
        (contfgl=contfgl,magicped=magicped)
    end
end

function pedinfo2magicped(pedinfo::Union{AbstractString,MateScheme};    
    popsize::Union{Nothing,Integer}=200,
    isfounderinbred::Bool=true,
    workdir::AbstractString = pwd())
    if isa(pedinfo,MateScheme)
        magicped = rand_magicped(pedinfo; popsize)        
    else
        if last(splitext(pedinfo))==".csv"
            magicped = readmagicped(pedinfo; workdir)
        else
            magicped = formmagicped(pedinfo,popsize)            
        end
        if isa(magicped.designinfo, Dict{String, DesignInfo}) 
            ped = MagicBase.magicped2ped(magicped; isfounderinbred)
            nf = ped.nfounder
            if ped.member[1:nf] == magicped.founderinfo[!,:individual]
                magicped.designinfo = ped
            else
                dict = Dict(ped.member[1:nf] .=> 1:nf)
                indices = [get(dict,i,missing) for i in magicped.founderinfo[!,:individual]]            
                any(ismissing.(indices)) && @error string("missing founders!")
                df = Pedigrees.ped2df(ped)
                df[1:nf,:] .= df[indices,:]
                magicped.designinfo = Pedigree(df)
            end
        end
    end
    magicped
end

function ispedfile(pedinfo)
    isa(pedinfo,AbstractString) && last(splitext(pedinfo)) == ".csv"
end

