
function magicimpute_founder!(magicgeno::MagicGeno;
    model::AbstractString="jointmodel",
	likeparameters::LikeParameters=LikeParameters(),   
	threshlikeparameters::ThreshLikeParameters=ThreshLikeParameters(),    
	priorlikeparameters::PriorLikeParameters=PriorLikeParameters(),    
	israndallele::Bool=true,
	isfounderinbred::Bool=true,				
	byfounder::Integer=0,	
	isrepeatimpute::Union{Nothing,Bool}=false, 
    nrepeatmin::Integer=3,
    nrepeatmax::Integer=6,     
	inputneighbor::Union{Nothing,AbstractDict}=nothing,
	isinferjunc::Union{Nothing, Bool} = nothing, 
	iscorrectfounder::Union{Nothing, Bool} = nothing,    
    isdelmarker::Bool = true,
    delsiglevel::Real = 0.01,    		
	skeletonsize::Union{Nothing,Integer} = nothing, 	
	isinfererror::Union{Nothing, Bool} = true,
	tukeyfence::Real=3.0,						
	minoutlier::Real=0.05, 
	isimputefounder::Union{Nothing,Bool}=nothing, 
	isallowmissing::Bool=true,
	isordermarker::Bool = !isnothing(inputneighbor),
	isspacemarker::Bool = !isnothing(inputneighbor) || isordermarker,
    trimcm::Real=20,
	trimfraction::Real=0.05,  #cM
	slidewin::Union{Nothing,Integer} = nothing, 
	slidewin_neighbor::Union{Nothing,Integer} = 200,
	orderactions::AbstractVector = ["inverse","permute"],  
    orderactions_neighbor::AbstractVector = ["inverse11","inverse01"],  
    inittemperature::Real= isordermarker ? 2.0 : 0.0,
    coolrate::Real=0.85,
    minaccept::Real=0.15,
	spacebyviterbi::Bool=false, 
	isparallel::Bool=true,	
	maxiter::Integer = 50,
    workdir::AbstractString=pwd(),
    tempdirectory::AbstractString = tempdir(),
	outstem::Union{Nothing,AbstractString}="outstem",
    outext::AbstractString=".vcf.gz",
    logfile::Union{Nothing,AbstractString,IO}= outstem*"_founder.log",
    verbose::Bool=true,
	more_verbose::Bool=false)
    starttime = time()
    io = MagicBase.set_logfile_begin(logfile, workdir, "magicimpute_founder!"; verbose)
	seqerror = MagicBase.get_seqerror(likeparameters)
	epsf = MagicBase.get_foundererror(likeparameters)
	epso = MagicBase.get_offspringerror(likeparameters)	
	MagicReconstruct.check_common_arg(model, epsf, epso, seqerror,
        workdir, tempdirectory)
	check_impute_arg(0.9,byfounder,
        delsiglevel,trimcm, trimfraction,minaccept,inittemperature, coolrate)
    model = MagicBase.reset_model(magicgeno.magicped,model;io,verbose)
    isparallel = isparallel && nprocs() > 1	&& (length(magicgeno.markermap) > 1 || nrepeatmax == nrepeatmin > 1)
    msg = string("list of options: \n",
        "model = ", model, "\n",
		"likeparameters = ", likeparameters, "\n",		
		"threshlikeparameters = ", threshlikeparameters, "\n",		
		"priorlikeparameters = ", priorlikeparameters, "\n",	
		"israndallele = ", israndallele,"\n",				
		"isfounderinbred = ", isfounderinbred,"\n",								
		"byfounder = ", byfounder, "\n",
		"isallowmissing = ", isallowmissing, "\n",
		"isrepeatimpute = ", isrepeatimpute, "\n",	
		"nrepeatmin = ", nrepeatmin, "\n",	
		"nrepeatmax = ", nrepeatmax, "\n",	
		"inputneighbor = ", isnothing(inputneighbor) ? "nothing" : string("dict of size ",length(inputneighbor)), "\n",
		"isinferjunc = ", isinferjunc, "\n",
		"iscorrectfounder = ", iscorrectfounder, "\n",		
        "isdelmarker = ", isdelmarker, "\n",
        "delsiglevel = ", delsiglevel, "\n",
        "isinfererror = ", isinfererror, "\n",				
		"tukeyfence = ", tukeyfence, "\n",								
		"minoutlier = ", minoutlier, "\n",								
		"isordermarker = ", isordermarker, "\n",
		"isspacemarker = ", isspacemarker, "\n",
        "trimcm = ", trimcm, "\n",
		"trimfraction = ", trimfraction, "\n",
		"skeletonsize = ", skeletonsize, "\n",
        "slidewin = ", slidewin, "\n",
		"slidewin_neighbor = ", slidewin_neighbor, "\n",
		"orderactions = ", orderactions, "\n",
		"orderactions_neighbor = ", orderactions_neighbor, "\n",
		"inittemperature = ", inittemperature, "\n",		
        "coolrate = ", coolrate, "\n",
		"minaccept = ", minaccept, "\n",
		"spacebyviterbi = ", spacebyviterbi, "\n",		
		"isparallel = ", isparallel, isparallel ? string("(nworkers=",nworkers(),")") : "", "\n",
		"maxiter = ",maxiter,"\n",
        "workdir = ",workdir,"\n",
        "tempdirectory = ",tempdirectory,"\n",
		"outstem = ",outstem,"\n",
		"outext = ",outext,"\n",
        "logfile = ",logfile,"\n",
        "verbose = ",verbose, "\n",
		"more_verbose = ", more_verbose)
    printconsole(io,verbose,msg)	
	# if isnothing(skeletonsize)
	# 	skeletonsize = model == "depmodel" ? 200 : 100
	# 	msg = string("reset skeletonsize=",skeletonsize)
	# 	printconsole(io,verbose,msg)
	# end
	if !isnothing(inputneighbor)
		msg = string("average #neighbors per marker =  ", round(mean(length.(values(inputneighbor))),digits=1))
		printconsole(io,verbose,msg)
	end
	if isnothing(isinferjunc)
		isinferjunc = false
		designinfo = magicgeno.magicped.designinfo
		if isa(magicgeno.magicped.designinfo,JuncDist)	
			if model in ["depmodel","indepmodel"]
				isinferjunc = isnothing(designinfo.mapexpansion) && isnothing(designinfo.j1122)
			else
				(; ibd, j1122, j1211, j1213, j1222,j1232) = designinfo
				isinferjunc = any(isnothing.([ibd, j1122, j1211, j1213, j1222,j1232]))
			end
		end
		msg = string("reset isinferjunc = ", isinferjunc)
		printconsole(io,verbose,msg)
	end
	if isinferjunc && isspacemarker
		isinferjunc = false
		msg = string("juncdist.mapexpansion and chrlen are indistinguishable. reset isinferjunc=false")
		@warn msg
		printconsole(io,false,string("Warning: ",msg))
	end
	if isinferjunc && in(model,["jointmodel","indepmodel"])
		isinferjunc = false
		msg = string("TODO: for jointmodel and indepmodel. reset isinferjunc=false")
		@warn msg
		printconsole(io,false,string("Warning: ",msg))
	end		
	MagicBase.reset_juncdist!(magicgeno.magicped,model;
		io,verbose,isfounderinbred,isinferjunc)				
	if isordermarker || !isfounderinbred
		nsnpmiss = 0
	else
		tused = @elapsed missgeno =  MagicBase.split_missprogeny!(magicgeno)				
		nsnpmiss = sum(size.(missgeno.markermap,1))		
		printconsole(io,verbose,msg)		
		if nsnpmiss == 0
			msg = "no markers with all offspring genotypes being missing"
			printconsole(io,verbose,msg)
		else
			msg = string("remove (and insert afterwards) ", nsnpmiss, " markers with all offspring genotypes being missing")			
			printconsole(io,verbose,msg)
			mem1 = round(Int, memoryuse()/10^6)
			missgenofile = getabsfile(workdir,string("temporary_", outstem,"_magicimpute_missgeno.csv.gz"))
			savemagicgeno(missgenofile,missgeno)
			missgeno = nothing
			GC.gc()
			mem2 = round(Int, memoryuse()/10^6)
			msg = string("split missgeno, tused=",round(tused,digits=1), 
				"seconds, mem=",mem1,"|",mem2,"MB")
			printconsole(io,verbose,msg)		
		end
	end	   
	MagicBase.rawgenoprob!(magicgeno; targets = ["founders"], seqerror, isfounderinbred)
    MagicBase.rawgenocall!(magicgeno; targets = ["founders"], callthreshold = 0.95, isfounderinbred)
    founderformat, offspringformat = MagicBase.setunphasedgeno!(magicgeno)		
	liketargetls = MagicBase.get_liketargetls(likeparameters)
	if all(occursin.(r"^GT",offspringformat))
        setdiff!(liketargetls,["seqerror","allelebalancemean","allelebalancedisperse","alleledropout"])
    end
    if model == "depmodel"
        setdiff!(liketargetls,["allelebalancemean","allelebalancedisperse","alleledropout"])
    end
	if isnothing(isinfererror) 		
		isinfererror = isspacemarker || in("AD", offspringformat) 
		printconsole(io,verbose, string("reset isinfererror=",isinfererror))
	end
	if isempty(liketargetls) && isinfererror
		msg = string("reset isinfererror = false for fixed likeparameters = ",likeparameters)
		printconsole(io, false, string("Warning: ",msg))
		@warn msg
		isinfererror = false
	end
	if isinfererror
		msg = string("inferring likelihood parameters ", liketargetls)
		printconsole(io,verbose,msg)
	end
	if isnothing(iscorrectfounder) 		
		iscorrectfounder = !in("AD", offspringformat) || (model == "depmodel") 
		printconsole(io,verbose, string("reset iscorrectfounder=",iscorrectfounder))
	# else
	# 	if iscorrectfounder && in("AD", offspringformat) && model != "depmodel" && !isinfererror 
	# 		msg = string("founder correction might not work well for sequence data without inferring allelic bias")
	# 		@warn msg
	# 		printconsole(io,false, "Warning: "*msg)
	# 	end
	end	
	if in("GT_phased",founderformat)
        if isfounderinbred
            @error string("founderformat=",founderformat, ", inconsistent with isfounderinbred = ",isfounderinbred)        
        end
    end
	if in("GT_haplo",offspringformat)
        @error string("offspringformat=",offspringformat, ", offspring have haplotypes, not genotypes")
    end
	if isnothing(isrepeatimpute)		
		if isimputefounder === false
			isrepeatimpute2 = false 
		else
			isrepeatimpute2 = get_isrepeat_impute(magicgeno; isfounderinbred,maxfmiss=0.5)		
		end
		printconsole(io,verbose, string("reset isrepeatimpute = ", isrepeatimpute2))	
	else
		isrepeatimpute2 = isrepeatimpute
		if isimputefounder === false && isrepeatimpute		
			isrepeatimpute2 = false 
			msg = string("reset isrepeatimpute = ", isrepeatimpute2, " since isimputefounder=",isimputefounder)
			@warn msg
			printconsole(io,false, "Warning: "*msg)		
		end
	end	
	MagicBase.info_magicgeno(magicgeno;io,verbose)		
	describe_phase_msg(io, verbose,offspringformat)
	# isimpute_afterrepeat = true
	if isrepeatimpute2		
		nrepeatimpute =  nrepeatmin == nrepeatmax  ? nrepeatmin : (nrepeatmin,nrepeatmax)
		# isspacemarker=false, isordermarker=false: assuming input input markermap is good enough and founder imputation is roobust to the inital genetic map
		# uickinfererror=true: skip estimations for some error rate parameters to speed up computation		
		partmagicgenols = MagicBase.splitby_connectedped(magicgeno);
		# TODO: 
		if length(partmagicgenols) == 1 || isordermarker 			
			printconsole(io,verbose, string("\nstart founder imputation with nrepeatmin = ", nrepeatmin, " and nrepeatmax = ",nrepeatmax))						
			magicimpute_founder_repeat!(magicgeno,nrepeatimpute;
				model, likeparameters, threshlikeparameters, priorlikeparameters, quickinfererror = false,
				israndallele, isfounderinbred, byfounder,  
				inputneighbor, isinferjunc, iscorrectfounder,isimputefounder, isallowmissing, 
				isdelmarker, isinfererror, isordermarker, isspacemarker, 
				delsiglevel, tukeyfence,  minoutlier, trimcm, trimfraction, skeletonsize, 
				slidewin,slidewin_neighbor, orderactions, orderactions_neighbor, 
				inittemperature, coolrate, minaccept, spacebyviterbi,
				isparallel, maxiter, workdir,tempdirectory,
				io, verbose,more_verbose
			)
		else			
			# length(partmagicgenols) > 1 && !isordermarker 			
			for partmagicgeno in partmagicgenols
				if isnothing(isrepeatimpute)
					isrepeatimpute_sub = get_isrepeat_impute(partmagicgeno; isfounderinbred,maxfmiss=0.5)		
					isrepeatimpute_sub || continue				
				end
				subpopls = unique(partmagicgeno.magicped.offspringinfo[!,:member])								
				printconsole(io,verbose, string("\nstart founder imputation with nrepeatmin = ", 
					nrepeatmin, " and nrepeatmax = ",nrepeatmax, " for subpopulations = ",subpopls))						
				MagicBase.info_magicgeno(partmagicgeno;io,verbose)		
				magicimpute_founder_repeat!(partmagicgeno,nrepeatimpute;
					model, likeparameters, threshlikeparameters, priorlikeparameters, quickinfererror=false, 
					israndallele, isfounderinbred, byfounder,  
					inputneighbor, isinferjunc, iscorrectfounder, isimputefounder, isallowmissing,
					isdelmarker, isinfererror, isordermarker, isspacemarker, 
					delsiglevel, tukeyfence, minoutlier, trimcm, trimfraction, skeletonsize, 
					slidewin,slidewin_neighbor, orderactions, orderactions_neighbor, 
					inittemperature, coolrate, minaccept,spacebyviterbi,                
					isparallel, maxiter, workdir,tempdirectory,
					io, verbose,more_verbose
				)
				# copy estimated foundergeno in each connected pouplation
				# warning: if isinfererror, estimated error rates are not copied
				# warning: deleted markers in subpopulations are kept 
				founders = magicgeno.magicped.founderinfo[!,:individual]
				founderdict = Dict(founders .=> eachindex(founders))
				for chr in eachindex(magicgeno.markermap)
					markers = magicgeno.markermap[chr][!,:marker]
					markerdict = Dict(markers .=> eachindex(markers))				
					snppos = [markerdict[i] for i in partmagicgeno.markermap[chr][!,:marker]]
					founderpos = [founderdict[i] for i in partmagicgeno.magicped.founderinfo[!,:individual]]					
					magicgeno.foundergeno[chr][snppos, founderpos] .= partmagicgeno.foundergeno[chr]					
					if !isfounderinbred
						magicgeno.markermap[chr][snppos,:founderformat] .= "GT_phased"
					end
				end			
			end	
			# delete markers with "GT_unphased" for outbred parents, or to consider mixed phased/unphased markers
			#  at a given marker, some parents' genotypes may be unphased! set them to phased genotypes if they are homozygoutes and missing if heterozygous			
			isfounderinbred || remove_GT_unphased!(magicgeno)
			# reduce memoryuse
			for partmagicgeno in partmagicgenols
				partmagicgeno.markermap = nothing
				partmagicgeno.foundergeno = nothing
				partmagicgeno.offspringgeno = nothing
			end
		end
	end		
	msg = string("start founder imputation without repeating")			
	MagicBase.reset_juncdist!(magicgeno.magicped,model;
		io,verbose,isfounderinbred,isinferjunc)		# reset_juncdist! is required
	isrepeatimpute2 && (msg = "\n"*msg)
	printconsole(io,verbose, msg)						
	nrepeatimpute = isrepeatimpute2 ? 1 : -1           	
	magicimpute_founder_repeat!(magicgeno,nrepeatimpute;
		model, likeparameters, threshlikeparameters, priorlikeparameters, quickinfererror=false, 
		israndallele, isfounderinbred, byfounder, 
		inputneighbor, isinferjunc, iscorrectfounder, isimputefounder, isallowmissing,	
		isdelmarker, isinfererror, isordermarker, isspacemarker, 
		delsiglevel, tukeyfence, minoutlier, trimcm, trimfraction, skeletonsize, 
		slidewin,slidewin_neighbor, orderactions, orderactions_neighbor, 
		inittemperature, coolrate, minaccept,spacebyviterbi,                 
		isparallel, maxiter, workdir,tempdirectory,
		io, verbose,more_verbose
	)	
	if nsnpmiss > 0							
		missgeno = readmagicgeno(missgenofile)
		# savemagicgeno("estgeno.csv.gz", magicgeno; workdir)		
		# savemagicgeno("missgeno.csv.gz", missgeno; workdir)
		# @info string("missgenofile=",missgenofile) 
		rm(missgenofile;force=true)		
		MagicBase.merge_missprogeny!(magicgeno,missgeno)
		msg = string("insert ", nsnpmiss, " markers with all offspring genotypes being missing")
		printconsole(io,verbose,msg)
		MagicBase.check_markerorder(magicgeno; io,verbose)
		if isinfererror
			for chr in 1:length(magicgeno.markermap)				
				for col in [:foundererror, :offspringerror,:seqerror,:allelebalancemean, :allelebalancedisperse,:alleledropout]
					ls = magicgeno.markermap[chr][!,col]
					b = ismissing.(ls)
					if !all(b)
						ls[b] .= mean(skipmissing(ls))					
						magicgeno.markermap[chr][!,col] .= float.(ls)
					end
				end
			end
		end
	end	
	if !isnothing(outstem)			
		# save peroffspringerr		
		offinfo = MagicBase.offspringinfo2df(magicgeno.magicped.offspringinfo)			
		offcols = occursin.(r"^peroffspringerror",names(offinfo))
		if any(offcols)
			nlargeerr = sum(Matrix(offinfo[!,offcols] .>  threshlikeparameters.peroffspringerror),dims=2)[:,1] 
			insertcols!(offinfo, size(offinfo,2)+1, "#LG_large" =>nlargeerr)
			islargeerr = nlargeerr .>= 1
			if any(islargeerr)				
				msg = string(sum(islargeerr)," offspring with too large peroffspringerror: ",offinfo[islargeerr,:individual])				
				@warn msg
				printconsole(io,false, "Warning: "*msg)                    
			else
				printconsole(io,verbose, "no offspring with too large peroffspringerror. ")
			end
			outputfile= string(outstem,"_peroffspringerror.csv")
			CSV.write(getabsfile(workdir,outputfile),offinfo)       
			msg = string("save peroffspringerror in ",outputfile)
			printconsole(io,verbose, msg)	 
		end		
		# save genodata
		outputfile= string(outstem,"_founder", outext)
		tused = @elapsed savegenodata(outputfile, magicgeno; workdir)
		msg = string("save founderimputed genofile: ",outputfile,
			", tused=",round(tused,digits=1),"s")
		printconsole(io,verbose,msg)
		# save founder deletion and correction
		for keystr in ["deletion","correction"]
			haskey(magicgeno.misc, keystr) || continue
			df = magicgeno.misc[keystr]
			if !isempty(df)
				outputfile= getabsfile(workdir, string(outstem,"_", keystr, ".csv"))
				CSV.write(outputfile, df)
				msg = string("save ", keystr, " file: ", relpath(outputfile,workdir))
				printconsole(io,verbose,msg)
			end
		end		
		# save mapfile and plot maps and errors
		if isordermarker || isspacemarker || isinfererror
			outmapfile = savemapfile(magicgeno; workdir, outstem)
			printconsole(io,verbose,string("save refinedmap file: ", outmapfile))    
			if isinfererror
				figerr = MagicBase.plotmarkererror(outmapfile;tukeyfence, workdir)
				if !isnothing(figerr)
					errfile = outstem*"_permarkererror.png"
					try 
						MagicBase.savefig(figerr,getabsfile(workdir,errfile))
					catch err
						msg = string("Could not savefig. ",err)
						@warn msg
						printconsole(io, verbose, "Warning: "*msg)
					end
				end
				if in("peroffspringerror",liketargetls)			
					peroff_errfile = string(outstem, "_peroffspringerror.csv") 
					figerr = MagicBase.plot_peroffspringerror(peroff_errfile; tukeyfence, workdir)
					if !isnothing(figerr)
						errfile = outstem*"_peroffspringerror.png"
						try 
							MagicBase.savefig(figerr,getabsfile(workdir,errfile))
						catch err
							msg = string("Could not savefig. ",err)
							@warn msg
							printconsole(io, verbose, "Warning: "*msg)
						end
					end
				end
			end
			if isordermarker || isspacemarker				
				markermap = reduce(vcat, magicgeno.markermap)
				figmap = plotmarkermap(markermap, copy(markermap); 					
					isphysmap = [true,false],
					cordigits = 2, 
					maplabels = ["Physical position(Mbp)", "Refined genetic position(cM)"],
				) 		
				if !isnothing(figmap)	
					figfile = outstem*"_compare_physmap.png"
					try 
						MagicBase.savefig(figmap,getabsfile(workdir,figfile))
					catch err
						msg = string("Could not savefig. ",err)
						@warn msg
						printconsole(io, verbose, "Warning: "*msg)
					end
				end
			end
        end				
	end
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, io, tused,"magicimpute_founder!"; verbose)
    magicgeno
end

function describe_phase_msg(io::Union{Nothing, IO},verbose::Bool,offspringformat::AbstractVector)
    msg = "keywords in print messages: \n"
    msg *= "\t#diff      ≡ #differences between proposal and current imputing\n"    
	msg *= "\tΔlogl      ≡ increase of log-likelihood due to proposal imputing\n"
	msg *= "\tstuck      ≡ measure of imputation being stuck\n"
	msg *= "\tbyhalf     ≡ if true, obtain proposal for one-half-chr and the other\n"
	msg *= "\t#correct_f ≡ #corrected founder genotypes (and set obsgeno missing)\n"
	msg *= "\t#mono      ≡ #monomorphic markers\n"
	msg *= "\t#del_mono  ≡ #deleted monomorphic markers\n"
	msg *= "\t#del_err   ≡ #deleted markers with too large error rates\n"
	msg *= "\t#del       ≡ #deleted markers via Vuong's test\n"
	msg *= "\tεo         ≡ offspring error rate per marker \n"
	msg *= "\tξo         ≡ offspring error rate per individual\n"
	if in("AD",offspringformat)	
		msg *= "\tεseq       ≡ sequence base error rate per marker\n"
		msg *= "\tABmean     ≡ sequence allelic bias per marker\n"
		msg *= "\tdisperse   ≡ sequence overdispersion per marker\n"
	end
	msg *= "\t#off_excl  ≡ #offspring temporarily excluded due to large ξo\n"
	msg *= "\tt          ≡ time used for each sub-step\n"	
	msg *= "\tmem        ≡ memory use before and after garbage collection"	
    printconsole(io, verbose,msg)
end

function remove_GT_unphased!(magicgeno::MagicGeno) 
	# only for outbred parents
	# delete markers with "GT_unphased" for outbred parents, or to consider mixed phased/unphased markers
	for chr in eachindex(magicgeno.markermap)
		chrmarkermap = magicgeno.markermap[chr]				
		b = chrmarkermap[!,:founderformat] .== "GT_phased"
		all(b) && continue
		keepat!(chrmarkermap,b)
		magicgeno.foundergeno[chr] = magicgeno.foundergeno[chr][b,:]
		magicgeno.offspringgeno[chr] = magicgeno.offspringgeno[chr][b,:]
	end
	# at a given marker, some parents' genotypes may be unphased! set them to phased genotypes if they are homozygoutes and missing if heterozygous			
	for chr in eachindex(magicgeno.markermap)
		chrmarkermap = magicgeno.markermap[chr]		
		all(chrmarkermap[!, :founderformat] .== "GT_phased") || @error string("outbred founders are not phased!")
		chrfgeno = magicgeno.foundergeno[chr]
		nf = size(chrfgeno,2)
		nsnp = size(chrmarkermap,1)
		for snp in 1:nsnp
			for j in 1:nf
				if isa(chrfgeno[snp,j],AbstractString) 
					gphased = split(chrfgeno[snp,j],"")					
					if allequal(gphased)
						chrfgeno[snp,j] = gphased
					else
						chrfgeno[snp,j] = ["N","N"]
					end								
				end
			end			
		end
	end	
	nothing
end

function get_isrepeat_impute_old(magicgeno::MagicGeno,chr::Integer;    
    magicprior::NamedTuple, 
    model::AbstractString,     
    isfounderinbred::Bool,
    byfounder::Integer=0,
    maxfmiss=0.5)
	fhaplosetpp = MagicImpute.calfhaploprior(magicgeno,chr)		
	fmissls = get_fmissls(fhaplosetpp) 	
    popmakeup = MagicReconstruct.calpopmakeup(magicgeno,chr,model, magicprior; isfounderinbred)    	
	founderformat = unique(magicgeno.markermap[chr][!,:founderformat])	
	isfounderphased = issubset(["GT_phased"], founderformat)
	defaultby = MagicBase.get_partition_defaultby(popmakeup, isfounderinbred,isfounderphased)
	findexlist = MagicBase.getfindexlist(byfounder,fmissls, popmakeup;defaultby,minfmiss=1e-4)	
	length(findexlist) > 1 && (!isfounderinbred || maximum(fmissls) > maxfmiss)
end


function get_isrepeat_impute(magicgeno::MagicGeno,chr::Integer;        
    isfounderinbred::Bool, maxfmiss=0.5)
	if isfounderinbred
		fhaplosetpp = MagicImpute.calfhaploprior(magicgeno,chr)		
		fmissls = get_fmissls(fhaplosetpp) 	    
		maximum(fmissls) > maxfmiss					
	else
		true
	end	
end

function get_isrepeat_impute(magicgeno::MagicGeno;     
    isfounderinbred::Bool,    
    maxfmiss=0.5)    
    nchr = length(magicgeno.markermap)
    any([get_isrepeat_impute(magicgeno,chr; isfounderinbred,maxfmiss) for chr in 1:nchr])
end



function get_magicgenofilemtx(magicgenofilels::AbstractVector,nrepeat::Integer)
	if length(magicgenofilels) <=1
		tsleepmax = 1
	else
		tsleepls = parse.(Float64,replace.(last.(split.(first.(MagicBase.split_allext.(basename.(magicgenofilels))),"_")),"tsleep"=>""))
		tsleep_step = mean(abs.(diff(tsleepls)))		
		tsleepmax = max(1,round(Int,maximum(tsleepls)+tsleep_step))
	end
    magicgenofilemtx = permutedims(reduce(hcat,[get_magicgenofilels(file,nrepeat; tsleepmax) for file in magicgenofilels]))
    magicgenofilemtx
end

function get_magicgenofilels(magicgenofile::AbstractString,nrepeat::Integer; tsleepmax::Real=0)
	filebase,fileext = MagicBase.split_allext(magicgenofile)
	filedir = dirname(filebase)
	filebase2 = basename(filebase)
	tsleepstr = last(split(filebase2,"_"))
	occursin("tsleep",tsleepstr) || @error string("unexpected filename=",filebase)
	filetsleep = parse(Float64,replace(tsleepstr,"tsleep"=>""))    
	newbase = joinpath(filedir,join(split(filebase2,"_")[1:end-1],"_"))
	ls = [string(newbase,"_tsleep", round(Int,filetsleep + (j-1)*tsleepmax),"_run",j,fileext) for j in 1:nrepeat]	
	ls
end


function magicimpute_founder_repeat!(magicgeno::MagicGeno,nrepeatimpute::Tuple;
    model::AbstractString="jointmodel",	
	inputneighbor::Union{Nothing,AbstractDict}=nothing,
    likeparameters::LikeParameters= isnothing(inputneighbor) ? LikeParameters(peroffspringerror=0.0) : LikeParameters(), 	
	threshlikeparameters::ThreshLikeParameters=ThreshLikeParameters(),    
	priorlikeparameters::PriorLikeParameters=PriorLikeParameters(),    
	quickinfererror::Bool=false,
	israndallele::Bool=true,
	isfounderinbred::Bool=true,				
	byfounder::Integer=0,		
	isinferjunc::Bool=false,
	iscorrectfounder::Bool=false,
	isimputefounder::Union{Nothing,Bool}=nothing, 
	isallowmissing::Bool,
    isdelmarker::Bool = true,    
	isinfererror::Bool = false,	
	isordermarker::Bool = !isnothing(inputneighbor),
	isspacemarker::Bool = !isnothing(inputneighbor) || isordermarker,
	delsiglevel::Real = 0.01,    		
	tukeyfence::Real=3.0,		
	minoutlier::Real=0.05, 					
    trimcm::Real=20,
	trimfraction::Real=0.05,  #cM
	skeletonsize::Union{Nothing,Integer} = nothing, 	
	slidewin::Union{Nothing,Integer} = nothing, 
	slidewin_neighbor::Union{Nothing,Integer} = 200,
	orderactions::AbstractVector = ["inverse","permute"],  
    orderactions_neighbor::AbstractVector = ["inverse11","inverse01"],  
    inittemperature::Real= isordermarker ? 2.0 : 0.0,
    coolrate::Real=0.85,
    minaccept::Real=0.15,
	spacebyviterbi::Bool=false, 
	isparallel::Bool=true,	
	maxiter::Integer = 50,
    workdir::AbstractString=pwd(),
    tempdirectory::AbstractString = tempdir(),	
    io::Union{Nothing,IO}=nothing, 
    verbose::Bool=true,
	more_verbose::Bool=false)
	if isinfererror && quickinfererror		
		(;foundererror,offspringerror,peroffspringerror,seqerror,allelebalancemean,allelebalancedisperse,alleledropout) = likeparameters
		if isnothing(peroffspringerror)
			peroffspringerror  = 0.0
			msg = string("peroffspringerror is reset to 0 for efficiently repeating founder imputation!")
			printconsole(io,verbose,msg)
		end
		# to add condition for AD data
		# if isnothing(allelebalancedisperse)
		# 	allelebalancedisperse  = 0.0
		# 	msg = string("allelebalancedisperse is reset to 0 for efficiently repeating founder imputation!")
		# 	printconsole(io,verbose,msg)
		# end
		# if isnothing(alleledropout)
		# 	alleledropout  = 0.0
		# 	msg = string("alleledropout is reset to 0 for efficiently repeating founder imputation!")
		# 	printconsole(io,verbose,msg)
		# end		
		likeparameters = LikeParameters(;foundererror,offspringerror,peroffspringerror,seqerror,allelebalancemean,allelebalancedisperse,alleledropout)
	end	
	magicprior=MagicReconstruct.calmagicprior(magicgeno,model; isfounderinbred)	
	#  split and save magicgeno by chromosomes
	# jldopen permissiond denied if tempdir() is in network drive    
	chridls = [i[1,:linkagegroup] for i in magicgeno.markermap]		
	jltempdir = mktempdir(tempdirectory; prefix="jl_magicimpute_founder_", cleanup=true)
	tempid = tempname(jltempdir,cleanup=false)		
	imputetempfilels = [string(tempid, "_impute_founder_",chrid, ".tmp") for chrid in chridls]
	nsnpls = [size(i,1) for i in magicgeno.markermap]
	chroo = MagicBase.get_chroo(nsnpls)	
	msg = string("chrid=>#snp: ",join([string(chridls[i],"=>",nsnpls[i]) for i in chroo],"|"))
	printconsole(io,verbose,msg)	
	jltempdir = mktempdir(tempdirectory; prefix="jl_magicimpute_founder_", cleanup=true)
	tempid = tempname(jltempdir,cleanup=false)	
	tused = @elapsed magicgenofilels = MagicBase.saveby_chromosome(magicgeno; 
		nworker=nworkers(), outstem = tempid)			
	mem1 = round(Int, memoryuse()/10^6)
	magicgeno.markermap = nothing
	magicgeno.foundergeno = nothing
	magicgeno.offspringgeno = nothing	
	GC.gc()
	GC.gc()
	mem2 = round(Int, memoryuse()/10^6)
	msg = string("saveby_chromosome, tused=",round(tused,digits=1), 
		"seconds, mem=",mem1,"|",mem2,"MB")
	printconsole(io,verbose,msg)	
	startt = time()	
    try
        if isparallel && nprocs()>1         
			tempid = tempname(jltempdir,cleanup=false)
			logfilels = [string(tempid, "_",chrid,"_run",run,".log") for chrid in chridls]
			try
		        pmap((x,y,z)->impute_refine_repeat_chr!(x,nrepeatimpute;
	                magicprior,
					model,israndallele, isfounderinbred,
					byfounder,
	                isinferjunc, iscorrectfounder,isimputefounder, isallowmissing,
					isdelmarker, delsiglevel,
					isspacemarker, trimcm, trimfraction,skeletonsize,
					isinfererror, likeparameters, threshlikeparameters, priorlikeparameters, tukeyfence, minoutlier, 															
					isordermarker, inputneighbor,slidewin, slidewin_neighbor,orderactions, orderactions_neighbor,
					inittemperature, coolrate, minaccept,spacebyviterbi, 
	                logio=y, maxiter, imputetempfile = z, 
					verbose,more_verbose),
	                magicgenofilels[chroo], logfilels[chroo],imputetempfilels[chroo])								
			catch err 				
				printconsole(io,verbose, string("ERROR: ",err))
				for (exc, bt) in current_exceptions()
					showerror(stdout, exc, bt)
				  	showerror(io, exc, bt)
					println(stdout)
				  	println(io)
				end
			finally
                if !isnothing(io)
					for file in logfilels
						if isfile(file)
							write(io,read(file))			            	
						end
			        end
				end
				rm.(logfilels,force=true)
			end
        else
			# chrls = sortperm([size(i,1) for i in magicgeno.markermap];rev=true)			
            for i in eachindex(magicgenofilels,imputetempfilels)				
				impute_refine_repeat_chr!(magicgenofilels[i],nrepeatimpute;
					magicprior,
					model,israndallele, isfounderinbred,
					byfounder,
	                isinferjunc, iscorrectfounder,isimputefounder, isallowmissing, 
					isdelmarker, delsiglevel,
					isspacemarker, trimcm, trimfraction,skeletonsize,
					isinfererror, likeparameters, threshlikeparameters, priorlikeparameters, tukeyfence, minoutlier, 															
					isordermarker, inputneighbor,slidewin,slidewin_neighbor,orderactions,orderactions_neighbor,
					inittemperature, coolrate, minaccept, spacebyviterbi,                   
					logio=io, maxiter, 
					imputetempfile = imputetempfilels[i], 
					verbose,more_verbose)				
            end
        end					
		msg = string("tused =", round(time()-startt,digits=1), "s for all ", length(magicgenofilels), " chromosomes ")
		MagicBase.printconsole(io,verbose,msg)		
		tused = @elapsed MagicBase.readby_chromosome!(magicgeno, magicgenofilels; workdir)
		mem1 = round(Int, memoryuse()/10^6)
		GC.gc()
		mem2 = round(Int, memoryuse()/10^6)		
		chrlen = round(sum([i[end,:poscm] - i[1,:poscm] for i in magicgeno.markermap]),digits=2)
		nsnp = sum(size.(magicgeno.markermap,1))
		msg = string("readby_chromosome! #markers=", nsnp, ", genome_length=", chrlen, "cM, tused=",round(tused,digits=1), 
			"seconds, mem=",mem1,"|",mem2,"MB")
		printconsole(io,verbose,msg)	
		nothing
    finally
		rm.(magicgenofilels,force=true)		        
		rm.(imputetempfilels,force=true)		
        rm(jltempdir,force=true,recursive=true)
    end	
end

function magicimpute_founder_repeat!(magicgeno::MagicGeno,nrepeatimpute::Integer;
    model::AbstractString="jointmodel",	
	inputneighbor::Union{Nothing,AbstractDict}=nothing,
    likeparameters::LikeParameters= isnothing(inputneighbor) ? LikeParameters(peroffspringerror=0.0) : LikeParameters(), 	
	threshlikeparameters::ThreshLikeParameters=ThreshLikeParameters(),    
	priorlikeparameters::PriorLikeParameters=PriorLikeParameters(),    
	quickinfererror::Bool=false,
	israndallele::Bool=true,
	isfounderinbred::Bool=true,				
	byfounder::Integer=0,	
	isinferjunc::Bool=false,
	iscorrectfounder::Bool=false,
	isimputefounder::Union{Nothing,Bool}=nothing, 
	isallowmissing::Bool, 
    isdelmarker::Bool = true,    
	isinfererror::Bool = false,	
	isordermarker::Bool = !isnothing(inputneighbor),
	isspacemarker::Bool = !isnothing(inputneighbor) || isordermarker,
	delsiglevel::Real = 0.01,    		
	tukeyfence::Real=3.0,			
	minoutlier::Real=0.05, 
    trimcm::Real=20,
	trimfraction::Real=0.05,  #cM
	skeletonsize::Union{Nothing,Integer} = nothing, 	
	slidewin::Union{Nothing,Integer} = nothing, 
	slidewin_neighbor::Union{Nothing,Integer} = 200,
	orderactions::AbstractVector = ["inverse","permute"],  
    orderactions_neighbor::AbstractVector = ["inverse11","inverse01"],  
    inittemperature::Real= isordermarker ? 2.0 : 0.0,
    coolrate::Real=0.85,
    minaccept::Real=0.15,
	spacebyviterbi::Bool=false, 
	isparallel::Bool=true,	
	maxiter::Integer = 50,
    workdir::AbstractString=pwd(),
    tempdirectory::AbstractString = tempdir(),	
    io::Union{Nothing,IO}=nothing, 
    verbose::Bool=true,
	more_verbose::Bool=false)
	if isinfererror && quickinfererror		
		(;foundererror,offspringerror,peroffspringerror,seqerror,allelebalancemean,allelebalancedisperse,alleledropout) = likeparameters
		if isnothing(peroffspringerror)
			peroffspringerror  = 0.0
			msg = string("peroffspringerror is reset to 0 for efficiently repeating founder imputation!")
			printconsole(io,verbose,msg)
		end
		# to add condition for AD data
		# if isnothing(allelebalancedisperse)
		# 	allelebalancedisperse  = 0.0
		# 	msg = string("allelebalancedisperse is reset to 0 for efficiently repeating founder imputation!")
		# 	printconsole(io,verbose,msg)
		# end
		# if isnothing(alleledropout)
		# 	alleledropout  = 0.0
		# 	msg = string("alleledropout is reset to 0 for efficiently repeating founder imputation!")
		# 	printconsole(io,verbose,msg)
		# end		
		likeparameters = LikeParameters(;foundererror,offspringerror,peroffspringerror,seqerror,allelebalancemean,allelebalancedisperse,alleledropout)
	end
	inputnrepeatimpute = nrepeatimpute
	nrepeatimpute = max(1, nrepeatimpute)
	magicprior=MagicReconstruct.calmagicprior(magicgeno,model; isfounderinbred)
	#  split and save magicgeno by chromosomes
	# jldopen permissiond denied if tempdir() is in network drive    
	chridls = [i[1,:linkagegroup] for i in magicgeno.markermap]
	jltempdir = mktempdir(tempdirectory; prefix="jl_magicimpute_founder_", cleanup=true)
	tempid = tempname(jltempdir,cleanup=false)		
	imputetempfilemtx = [begin 
		runstr=nrepeatimpute > 1 ? string("_run",run) : "_run"; 
		string(tempid, "_impute_founder_",chrid, runstr, ".tmp") 
	end for chrid in chridls, run in 1:nrepeatimpute]			
	nsnpls = [size(i,1) for i in magicgeno.markermap]
	chroo = MagicBase.get_chroo(nsnpls)	
	msg = string("chrid=>#snp: ",join([string(chridls[i],"=>",nsnpls[i]) for i in chroo],"|"))
	printconsole(io,verbose,msg)	
	tempid = tempname(jltempdir,cleanup=false)	
	tused = @elapsed magicgenofilels = MagicBase.saveby_chromosome(magicgeno; 
		nworker=nworkers(), outstem = tempid)	
	magicgenofilemtx = get_magicgenofilemtx(magicgenofilels,nrepeatimpute)	
	# copy magicgeno in magicgenofilels, that is, the first column of magicgenofilemtx
	for i in 1:length(magicgenofilels)
		for j in 1:nrepeatimpute
			cp(magicgenofilels[i],magicgenofilemtx[i,j]; force=true)
		end
	end	
	mem1 = round(Int, memoryuse()/10^6)
	magicgeno.markermap = nothing
	magicgeno.foundergeno = nothing
	magicgeno.offspringgeno = nothing	
	GC.gc()
	GC.gc()
	mem2 = round(Int, memoryuse()/10^6)
	msg = string("saveby_chromosome, tused=",round(tused,digits=1), 
		"seconds, mem=",mem1,"|",mem2,"MB")
	printconsole(io,verbose,msg)	
	startt = time()
	loglikemtx = zeros(size(magicgenofilemtx)...)	
	if inputnrepeatimpute == -1
		repeatrunmtx = [nothing for _ in 1:length(chridls), _ in 1:nrepeatimpute]
	else
		repeatrunmtx = [j for _ in 1:length(chridls), j=1:nrepeatimpute]
	end
    try
        if isparallel && nprocs()>1         
			tempid = tempname(jltempdir,cleanup=false)
			logfilemtx = [string(tempid, "_",chrid,"_run",run,".log") for chrid in chridls,run in 1:nrepeatimpute]							
			try
		        res = pmap((x,y,z,w)->impute_refine_chr!(x;
	                magicprior,
					model,israndallele, isfounderinbred,
					byfounder,
	                isinferjunc, iscorrectfounder,isimputefounder, isallowmissing,
					isdelmarker, delsiglevel,
					isspacemarker, trimcm, trimfraction,skeletonsize,
					isinfererror, likeparameters, threshlikeparameters, priorlikeparameters, tukeyfence, minoutlier, 															
					isordermarker, inputneighbor,slidewin, slidewin_neighbor,orderactions, orderactions_neighbor,
					inittemperature, coolrate, minaccept,spacebyviterbi, 
	                imputetempfile=y,logio=z, repeatrun=w, maxiter, 
					verbose,more_verbose),
	                magicgenofilemtx[chroo,:],imputetempfilemtx[chroo,:], logfilemtx[chroo,:],repeatrunmtx[chroo,:])				
				loglikemtx[chroo,:] .= sum.(first.(res))				
			catch err 				
				printconsole(io,verbose, string("ERROR: ",err))				
				for (exc, bt) in current_exceptions()
					showerror(stdout, exc, bt)
				  	showerror(io, exc, bt)
					println(stdout)
				  	println(io)
				end
			finally
                if !isnothing(io)
					for file in permutedims(logfilemtx)
						if isfile(file)
							write(io,read(file))			            	
						end
			        end
				end
				rm.(logfilemtx,force=true)
			end
        else
			# chrls = sortperm([size(i,1) for i in magicgeno.markermap];rev=true)			
            for i in eachindex(magicgenofilemtx)				
				chrres = impute_refine_chr!(magicgenofilemtx[i];
					magicprior,
					model,israndallele, isfounderinbred,
					byfounder,
	                isinferjunc, iscorrectfounder,isimputefounder, isallowmissing, 
					isdelmarker, delsiglevel,
					isspacemarker, trimcm, trimfraction,skeletonsize,
					isinfererror, likeparameters, threshlikeparameters, priorlikeparameters, tukeyfence, minoutlier, 															
					isordermarker, inputneighbor,slidewin,slidewin_neighbor,orderactions,orderactions_neighbor,
					inittemperature, coolrate, minaccept, spacebyviterbi, 
                    imputetempfile=imputetempfilemtx[i],repeatrun = repeatrunmtx[i], 					
					logio=io, maxiter, 
					verbose,more_verbose)				
				loglikemtx[i] = sum(chrres[1])
            end
        end			
		if inputnrepeatimpute >= 1	
			msg = "best runs: \n"
			for i in eachindex(chridls)
				loglikels = loglikemtx[i,:]
				bestlike,bestrun = findmax(loglikels)
				magicgenofilels[i] = magicgenofilemtx[i,bestrun]				
				bestgeno = readmagicgeno(magicgenofilemtx[i,bestrun])				
				bestfhaplo = get_chr_fhaplodf(bestgeno,1)
				founders = 1:size(bestgeno.magicped.founderinfo,1)
				fgenodiff = zeros(nrepeatimpute)
				fphasediff = zeros(nrepeatimpute)
				for j in 1:nrepeatimpute
					if j != bestrun
						rungeno = readmagicgeno(magicgenofilemtx[i,j])				
						runfhaplo = get_chr_fhaplodf(rungeno,1)
						fgenodiff[j] = caldiff_fhaplo(bestfhaplo,runfhaplo,founders; isfounderinbred)					
						if !isfounderinbred
							fphasediff[j] = calphasediff_fhaplo(bestfhaplo,runfhaplo,founders)							 							
						end
					end
				end
				loglikels[bestrun] = -Inf
				nbrlike,nbrrun = findmax(loglikels)
				if fgenodiff[nbrrun] >  0					
					nbrgeno = readmagicgeno(magicgenofilemtx[i,nbrrun])				
					nbrfhaplo = get_chr_fhaplodf(nbrgeno,1)
					isdiff = isdiff_fhaplo(bestfhaplo,nbrfhaplo,founders; isfounderinbred)							
					bestgeno.markermap[1][!,:marker] == bestfhaplo[!,:marker] || @error "inconsisent markers"					
					ndiff = sum(isdiff)					
				else
					ndiff = 0
				end
				msg *= string("chr=",chridls[i],", bestrun=",bestrun, ", bestlogl=",  round(bestlike,digits=1), 
					", lolhis = ",round.(loglikemtx[i,:],digits=1), ", genodiff=",round.(fgenodiff,digits=3))
				if !isfounderinbred
					msg *= string(", phasediff=", round.(fphasediff,digits=3))
				end
				if fgenodiff[nbrrun] >  0
					msg *= string(", #diff=",ndiff, " between bestrun ",bestrun, " and neighborrun ",nbrrun)
				end	
				if i < length(chridls)
					msg *= string("\n")
				end
			end
			MagicBase.printconsole(io,verbose,msg)		
		else 
			# inputnrepeatimpute == -1
			nrepeatimpute == 1 || @error string("unexepcted nrepeatimpute=",nrepeatimpute)
			magicgenofilels .= magicgenofilemtx[:,1]
		end		
		msg = string("tused =", round(time()-startt,digits=1), "s for all ", length(magicgenofilels), " chromosomes ")
		MagicBase.printconsole(io,verbose,msg)		
		tused = @elapsed MagicBase.readby_chromosome!(magicgeno, magicgenofilels; workdir)
		mem1 = round(Int, memoryuse()/10^6)
		GC.gc()
		mem2 = round(Int, memoryuse()/10^6)		
		chrlen = round(sum([i[end,:poscm] - i[1,:poscm] for i in magicgeno.markermap]),digits=2)
		nsnp = sum(size.(magicgeno.markermap,1))
		msg = string("readby_chromosome! #markers=", nsnp, ", genome_length=", chrlen, "cM, tused=",round(tused,digits=1), 
			"seconds, mem=",mem1,"|",mem2,"MB")
		printconsole(io,verbose,msg)	
		nothing
    finally
		rm.(magicgenofilemtx,force=true)		
        rm.(imputetempfilemtx,force=true)		
        rm(jltempdir,force=true,recursive=true)
    end	
end



function savemapfile(magicgeno::MagicGeno; 
	workdir::AbstractString=pwd(),    
	outstem::Union{Nothing,AbstractString}="outstem")
	markermap = reduce(vcat, magicgeno.markermap)
	markermap[!,:poscm] .= round.(markermap[!,:poscm],digits=5)
	outmapfile = string(outstem, "_map.csv")            
	outmapfile2 = getabsfile(workdir,outmapfile)
	open(outmapfile2,"w") do mapio
		descripls = [
			"marker, marker ID", 
			"linkagegroup, linkage group ID", 
			"poscm, marker position in centiMorgan",
			"physchrom, physical map chromosome ID for the marker",
			"physposbp, physical map position (in base pair) for the marker",
			"founderformat,format of founder genotypes at the marker",
			"offspringformat, format of offspring genotypes at the marker",
			"foundererror, allelic error rate in founders at the marker",
			"offsprignerror, allelic error rate in offspring at the marker",
			"seqeerror, sequence based error rate at the marker",
			"allelebalancemean, mean sequence allelic balance among offspring at the marker",
			"allelebalancedisperse, overdispersion for the distribution of sequence allelic balance among offspring at the marker",
			"alleledropout, 0- and 1-inflation for the distribution of sequence allelic balance among offspring at the maker"
		]
		msg = ""
		for i in eachindex(descripls)
			msg *= string("##col_", i,",", descripls[i],"\n")
		end
		write(mapio, msg)				
		CSV.write(mapio, markermap; delim=',', missingstring="NA",header=true,append=true)        
	end
	outmapfile
end

