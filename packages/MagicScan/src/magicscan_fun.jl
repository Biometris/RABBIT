"""
    magicscan(ancestryfile, phenofile; commentstring="##",kwargs...)

perform genomics scan of QTL for multiparental populations.
See `magicreconstruct` for generating ancestryfile by
haplotype reconstruction.

# Positional arguments

`ancestryfile::AbstractString` specifies ancestry file resulting `magicreconstruct`

`phenofile::AbstractString` specifies phenotypic file.

# Keyword arguments

`equation::Union{Nothing,FormulaTerm} = nothing` speficies linear model equation. By defulat,
StatsModel.@formula(y ~ 1), where y is the last column name in phenofile.

`thresholds::Union{Nothing,AbstractVector} = nothing`: list of thresholds used in plotting
scanning profile.

`nperm::Integer=200`: number of permutations of phenotypes  for calculating thresholds
that are not specified.

`siglevels::AbstractVector = [0.05]`: significance levels for calculating thresholds by permutations.

`islog10p::Bool=true`: if true, the profile refers to -log10P, and LOD otherwise.

`missingstring=["NA","missing"]`: string denotes a missing phenotypic value.

`commentstring::AbstractString` specifies the lines beginning with commentstring are ignored in genofile
or pedfile given by `pedinfo`.

`workdir::AbstractString=pwd()`: working directory for input and output files.

`outstem::Union{Nothing,AbstractString}="outstem"` specifies the stem of filename saving magicancestry.
See [`MagicBase.savemagicancestry`](@ref) for the description of outputfile "outstem_magicancestry.csv.gz".

`logfile::Union{Nothing,AbstractString,IO}= (isnothing(outstem) ? nothing : string(outstem,"_magicreconstruct.log"))`:
log file or IO for writing log. If it is nothing, no log file.

`verbose::Bool=true`: if true, print details on the stdout.

# Examples
```julia-repl
julia> magicscan("magicancestry.csv.gz","pheno.csv")
```
"""
function magicscan(ancestryfile::AbstractString,phenofile::AbstractString;
    equation::Union{Nothing,FormulaTerm} = nothing,
    thincm::Real = 0.0, 
    thresholds::Union{Nothing,AbstractVector} = nothing,
    nperm::Integer=200, # if thresholds is nothing, calculate thresholds by permuting phenotes.
    siglevels::AbstractVector = [0.05], # if thresholds is nothing, calculate thresholds at siglevels.
    islog10p::Bool=true, # if islog10p = true, profile refers to -log10P, and LOD otherwise.    
    missingstring=["NA","missing"],
    commentstring::AbstractString="##",
    workdir::AbstractString = pwd(),
    outstem::AbstractString="outstem",
    logfile::Union{AbstractString,IO}= string(outstem,"_magicscan.log"),
    verbose::Bool=true)
    starttime0 = time()
     if typeof(logfile) <: AbstractString
        io=open(MagicBase.getabsfile(workdir,logfile), "w+")
    else
        io=logfile
    end
    msg=string("Julia version: ",VERSION)
    MagicBase.printconsole(io,verbose,msg)
    MagicBase.printpkgst(io,verbose,"MagicScan")
    MagicBase.printconsole(io,verbose,string(repeat("=",33),"magicscan",repeat("=",34)))
    msg = string("magicscan, logfile=", logfile,
        ", ", round(Dates.now(),Dates.Second))
    MagicBase.printconsole(io,verbose,msg)
    MagicBase.printconsole(io,verbose,string("nthreads=",Threads.nthreads()))    
    msg = string("list of args: \n",
            "ancestryfile = ", ancestryfile, "\n",
            "phenofile = ", phenofile, "\n",
            "equation = ", equation, "\n",
            "thresholds = ", thresholds, "\n",
            "nperm = ", nperm, "\n",
            "siglevels = ", siglevels, "\n",
            "islog10p = ", islog10p, "\n",
            "missingstring = ", missingstring, "\n",
            "commentstring = ", commentstring, "\n",
            "outstem = ", outstem, "\n",
            "workdir = ", workdir, "\n",
            "logfile = ", logfile, "\n",
            "verbose = ", verbose)
    MagicBase.printconsole(io,verbose,msg)    
    outstem *= "_magicscan"
    # reading phenodf
    starttime = time()
    phenofile2 = getabsfile(workdir, phenofile)
    if !isfile(phenofile2) 
        msg = string(phenofile2, " does not exist!")
        MagicBase.printconsole(io,verbose,"Error: "*msg)
        error(msg)
    end
    #  missingstring is modified! after CSV.read
    phenodf = CSV.read(phenofile2,DataFrame;
        comment=commentstring, missingstring = copy(missingstring))
    if isnothing(equation)
        y = Symbol(last(names(phenodf)))
        f = Term(y) ~ ConstantTerm(1)
        MagicBase.printconsole(io,verbose,string("equation = ",f))
    else
        f = equation
    end    
    phenoname = string(f.lhs)
    filter!(row -> !ismissing(row[phenoname]),phenodf)
    phenodf[!,phenoname] .*= 1.0
    f = apply_schema(f, schema(f, phenodf))
    yresp,xpred = modelcols(f,phenodf)    
    # describe(io,yresp)
    # describe(yresp)
    try 
        his_pheno = histogram(yresp;legend=false,
            xlabel=string(f.lhs), 
            ylabel="frequency",
            left_margin=20Plots.mm,
            bottom_margin=15Plots.mm,                
        )
        outfile = string(outstem,"_response.png")
        outfile2 = getabsfile(workdir,outfile)
        savefig(his_pheno ,outfile2)
    catch
        @warn string("Could not plot phenotype histogram")
    end 
    # reading haploprob
    ancestryfile2 = getabsfile(workdir, ancestryfile)
    if !isfile(ancestryfile2) 
        msg = string(ancestryfile2, " does not exist!")
        MagicBase.printconsole(io,verbose,"Error: "*msg)
        error(msg)
    end
    magicancestry = readmagicancestry(ancestryfile;workdir)
    # check consistency of offspring and get haploprob
    phenooffspring = phenodf[!,1]
    genooffspring = magicancestry.magicped.offspringinfo[!,:individual]
    commonoffspring = intersect(phenooffspring,genooffspring)
    if isempty(commonoffspring)
        msg = "no common offspring between ancestryfile and phenofile"
        MagicBase.printconsole(io,false,string("Error: ",msg))
        error(msg)
    end
    if commonoffspring != phenooffspring
        d = setdiff(phenooffspring,genooffspring)
        if isempty(d)
            msg = string("all ", length(phenooffspring), " offspring in phenofile are also in ancestryfile")            
            MagicBase.printconsole(io,verbose,msg)
        else
            msg = string("ignore ", length(d), " out of ", length(phenooffspring), " offspring in phenofile but not in ancestryfile")
            MagicBase.printconsole(io,false,string("Warning: ",msg))
            @warn msg
        end
        dict = Dict(phenooffspring .=> 1:length(phenooffspring))
        offls = [dict[i] for i in commonoffspring]
        yresp = yresp[offls]
        xpred = xpred[offls,:]
    end
    if commonoffspring == genooffspring
        offls = 1:length(commonoffspring)
    else
        d = setdiff(genooffspring, phenooffspring)
        if isempty(d)
            msg = string("all ", length(genooffspring), " offspring in ancestryfile are also in phenofile")            
            MagicBase.printconsole(io,verbose,msg)            
        else
            msg = string("ignore ", length(d), " out of ", length(genooffspring), " offspring in ancestryfile but not in phenofile")
            MagicBase.printconsole(io,false,string("Warning: ",msg))
            @warn msg
        end
        dict = Dict(genooffspring .=> 1:length(genooffspring))
        offls = [dict[i] for i in commonoffspring]
    end
    nmarker_bef = sum(size.(magicancestry.markermap,1))    
    MagicBase.thinmagicancestry!(magicancestry; thincm) 
    nmarker_aft = sum(size.(magicancestry.markermap,1))   
    msg = string("#offspring=",size(phenodf,1),
        ", #marker_bef=",nmarker_bef,
        ", #marker_aft=",nmarker_aft,
        ", thincm=",thincm,"cM")
    MagicBase.printconsole(io,verbose,msg)
    haploprob = get_haploprob(magicancestry, offls)    
    msg = string("tused in reading data: ",
        round(time()-starttime,digits=2), "s")
    MagicBase.printconsole(io,verbose,msg)
    # scanning
    starttime = time()    
    loglikels,lodls,log10pls = scan_lm(yresp,xpred,haploprob);
    profile = reduce(vcat, magicancestry.markermap)[!,1:5]
    insertcols!(profile, size(profile,2)+1,
        :loglike =>reduce(vcat,loglikels),
        :lod =>reduce(vcat,lodls),
        :neg_log10p =>reduce(vcat,log10pls))
    outfile =string(outstem,"_profile.csv")
    outfile2 = getabsfile(workdir,outfile)    
    open(outfile2,"w") do outio
        descripls = ["marker, marker ID", 
            "linkagegroup, linkage group ID", 
            "poscm, marker position in centi-Morgan",
            "physchrom, physical map chromosome ID",
            "physposbp, physical map postion in base pair",
            "loglike, log likiehood for phenotypic data at the marker",
            "lod, LOD score for the liklihood ratio test at the marker",
            "neg_log10p, minus 10-based log of the P-value for the likelihood ratio test at the marker"
        ]
        for i in eachindex(descripls)
            write(outio, string("##col_",i, ", ", descripls[i], "\n"))
        end
        missingstring2 = isa(missingstring,AbstractVector) ? first(missingstring) : missingstring        
        CSV.write(outio, profile; delim=",", header=true, append=true, missingstring = missingstring2)
    end
    msg = string("tused in scanning QTL: ",
        round(time()-starttime,digits=3), "s")
    MagicBase.printconsole(io,verbose,msg)
    # calcuculating threshold
    starttime = time()
    outfile =string(outstem,"_perm.csv")
    outfile2 = getabsfile(workdir,outfile)
    if isnothing(thresholds)
        subpop2off = MagicBase.get_subpop2offspring(magicancestry.magicped;isindex=false)
        dict = Dict(commonoffspring .=> 1:length(commonoffspring))
        pop2off = Dict([begin
            intersect!(offs,commonoffspring)
            pop => [dict[i] for i in  offs]
        end for (pop,offs) in subpop2off])
        thesholddf = permute_threshold(yresp,xpred,haploprob,pop2off; nperm)
        col = islog10p ? :neg_log10p : :lod
        thresholds = [round(quantile(thesholddf[!,col],1-siglevel),digits=4)
         for siglevel in siglevels]
        msg = string("tused in calculating thresholds: ",
            round(time()-starttime,digits=2), "s\n",
            "\t", col, "-thresholds = ",thresholds, " at ", siglevels)
        MagicBase.printconsole(io,verbose,msg)
        open(outfile2,"w") do thresh_io
            msg = string("##", col,"-thresholds = ",thresholds,
                " at ", siglevels)
            write(thresh_io,msg,"\n")
            descripls = ["permute, index of permutation", 
                "loglike, log likelihood for th peak of scan profile obtained from the permuted phenotypes", 
                "lod, LOD score for th peak of scan profile obtained from the permuted phenotypes", 
                "neg_log10p, minus 10-based log of the P-value for the peak of scan profile obtained from the permuted phenotypes", 
            ]
            for i in eachindex(descripls)
                write(thresh_io, string("##col_",i, ", ", descripls[i], "\n"))
            end        
            CSV.write(thresh_io, thesholddf; header=true,append=true)
        end
        starttime = time()
    else
        isfile(outfile2) && rm(outfile2)
    end
    # save and plot
    try 
        manhattan = plot_manhattan(profile;thresholds,islog10p)
        qq = plot_qqpval(profile)
        manhattan_qq = plot(manhattan,qq;
            layout=(2,1),
            size = (1500, 800), 
            left_margin=20Plots.mm,
            bottom_margin=15Plots.mm,                
        )
        outfile = string(outstem,"_manhattan_qq.png")
        outfile2 = getabsfile(workdir,outfile)
        savefig(manhattan_qq ,outfile2)
    catch err
        @warn string("Could not plot manhattan. ", err)
    end 
    peak = getpeak(profile;CI_probs = 1.0 .- siglevels)
    peak13= peak[1:3,2]
    peak133 = tryparse(Float64,peak13[3])
    if !isnothing(peak133)
        peak13[3] = string(round(peak133,digits=4))
    end
    msg = string("peak=",join(peak13,"|"),", lod=", peak[4,2], ", -log10P=", peak[5,2])
    MagicBase.printconsole(io,verbose,msg)
    outfile = string(outstem,"_peak.csv")
    outfile2 = getabsfile(workdir,outfile)
    open(outfile2,"w") do outio
        descripls = ["id, variable ID && CI0.95 = 0.95 confidence interval for the peak && CI0.95_marker = CI in terms of nearest marker ID", 
            "peak, peak of scan profile"             
        ]
        for i in eachindex(descripls)
            write(outio, string("##col_",i, ", ", descripls[i], "\n"))
        end        
        CSV.write(outio, peak; delim=",", header=true, append=true)
    end
    msg = string("tused in saving and plotting: ",
        round(time()-starttime,digits=2), "s")
    MagicBase.printconsole(io,verbose,msg)
    MagicBase.printconsole(io,verbose,string("End, ", round(Dates.now(),Dates.Second),", tused = ",
        round(time()-starttime0,digits=2), " seconds by magicscan"))
    MagicBase.printconsole(io,verbose,repeat("=",76))
    if typeof(logfile) <: AbstractString
        close(io)
    elseif typeof(logfile) <: IO
        flush(io)
    end
    peak
end
