
function magicimpute_offspring!(magicgeno::MagicGeno;
	model::AbstractString="jointmodel",
	likeparam::LikeParam=LikeParam(),   
	israndallele::Bool=true,
	isfounderinbred::Bool=true,		
	threshimpute::Real=0.7,			
	phasealg::AbstractString="unphase", 
	isparallel::Bool=true,
    workdir::AbstractString=pwd(),
	tempdirectory::AbstractString = tempdir(),
	outstem::Union{Nothing,AbstractString}="outstem",
	outext::AbstractString=".vcf.gz",
    logfile::Union{Nothing,AbstractString,IO}= outstem*"_offspring.log",
    verbose::Bool=true)
	starttime = time()
	io = MagicBase.set_logfile_begin(logfile, workdir, "magicimpute_offspring!"; verbose)
	baseerror = MagicBase.get_likeproperty(likeparam, :baseerror)
	epsf = MagicBase.get_likeproperty(likeparam, :foundererror)
	epso = MagicBase.get_likeproperty(likeparam, :offspringerror)		
	MagicReconstruct.check_common_arg(model,epsf, epso, baseerror,
        workdir, tempdirectory)
	if !(0<=threshimpute<=1)
        @error string("threshimpute=", threshimpute, " is not in [0,1]")
    end
	if !in(phasealg, ["viterbi","forwardbackward"]) && !occursin(r"^unphas",phasealg)
		@error string("unknown phasealg=",phasealg)
	end
	model = MagicBase.reset_model(magicgeno.magicped,model;io,verbose)    	
	isparallel = isparallel && nprocs() > 1	&& length(magicgeno.markermap) > 1
	msg = string("list of options: \n",
        "model = ", model, "\n",
		"likeparam = ", likeparam, "\n",		
		"isfounderinbred = ", isfounderinbred,"\n",		
		"threshimpute = ", threshimpute,"\n",		
		"phasealg = ", phasealg,"\n",
        "isparallel = ", isparallel, isparallel ? string("(nworkers=",nworkers(),")") : "", "\n",
        "workdir = ",workdir,"\n",
		"tempdirectory = ",tempdirectory,"\n",
		"outstem = ",outstem,"\n",
		"outext = ",outext,"\n",
        "logfile = ",logfile,"\n",
        "verbose = ",verbose)
    printconsole(io,verbose,msg)
	MagicBase.reset_juncdist!(magicgeno.magicped,model;io,verbose,isfounderinbred)			
	MagicBase.rawgenoprob!(magicgeno; targets = ["founders"], baseerror, isfounderinbred)
    MagicBase.rawgenocall!(magicgeno; targets = ["founders"], callthreshold = 0.95, isfounderinbred)
	founderformat, offspringformat = MagicBase.setunphasedgeno!(magicgeno)	
	if in("GT_unphased",founderformat)
		@error string("founder genotypes are unphased! use MagicImpute to phase parents")
	end
	if in("GT_phased",founderformat)
        if isfounderinbred
            @error string("founderformat=",founderformat, ", inconsistent with isfounderinbred = ",isfounderinbred)        
        end
    end
	if in("GT_haplo",offspringformat)
        @error string("offspringformat=",offspringformat, ", offspring have haplotypes, not genotypes")
    end    
	MagicBase.info_magicgeno(magicgeno;io,verbose)	
	magicprior=MagicReconstruct.calmagicprior(magicgeno,model; isfounderinbred)	
	nsnpls = [size(i,1) for i in magicgeno.markermap]
	chroo = MagicBase.get_chroo(nsnpls)	
	#  split and save magicgeno by chromosomes
	# jldopen permissiond denied if tempdir() is in network drive    	
	nchr = length(magicgeno.markermap)
	jltempdir = mktempdir(tempdirectory; prefix="jl_magicimpute_offspring_", cleanup=true)
	tempid = tempname(jltempdir,cleanup=true)
    tempfilels = [string(tempid, "_impute_offspring_chr",chr,".tmp") for chr in 1:nchr]		
	tempid = tempname(jltempdir,cleanup=true)
	tused = @elapsed genofilels = MagicBase.saveby_chromosome(magicgeno; nworker = nworkers(),
		workdir=jltempdir, outstem=tempid)
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
			tempid = tempname(jltempdir,cleanup=true)
			logfilels = [string(tempid, "_chr",chr,".log") for chr in 1:nchr]				
			try
		        pmap((x,y,z)->impute_offspring_chr!(x, model,magicprior; likeparam, threshimpute,
					israndallele, isfounderinbred, phasealg, tempfile=y, logio=z,verbose),
					genofilels[chroo],tempfilels[chroo],logfilels[chroo])
			catch err
				msg = string("ERROR: ", err)
				printconsole(io,verbose,msg)
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
			            	flush(io)
						end
					end
				end
				rm.(logfilels,force=true)
				@everywhere GC.gc() # force garbage collection in all workers			
				@everywhere GC.gc() 
			end					
		else
			for chr in eachindex(genofilels)
				impute_offspring_chr!(genofilels[chr], model,magicprior; likeparam, threshimpute, 
					israndallele, isfounderinbred,phasealg, tempfile = tempfilels[chr], logio=io,verbose)
			end
		end
		msg = string("tused =", round(time()-startt,digits=1), "s for all chromosomes")
		MagicBase.printconsole(io,verbose,msg)	
		tused = @elapsed MagicBase.readby_chromosome!(magicgeno, genofilels; workdir)
		mem1 = round(Int, memoryuse()/10^6)
		GC.gc()
		mem2 = round(Int, memoryuse()/10^6)		
		msg = string("readby_chromosome!, tused=",round(tused,digits=1), 
			"seconds, mem=",mem1,"|",mem2,"MB")
		printconsole(io,verbose,msg)
		if !isnothing(outstem)
			if model != "depmodel"
				try 
					genofreqfile, figfile = MagicBase.plotgenofreq(magicgeno; 
						boundaryline = (1.5, :gray,:dash),
    					plotsize=(1500,1200),
    					minmaf_subpop=0.01, outstem,workdir)										
					MagicBase.printconsole(io,verbose,string("save in ", 
						isempty(dirname(outstem)) ? relpath(genofreqfile,workdir) : genofreqfile))
					MagicBase.printconsole(io,verbose,string("save in ", 
						isempty(dirname(outstem)) ? relpath(figfile,workdir) : figfile))
				catch err
					@warn string(err, ". Could not plot genofreq")
				end
			end
			outfile= string(outstem,"_geno", outext)
			outfile2 = getabsfile(workdir,outfile)
			tused = @elapsed save_imputed_gtgp(outfile2, magicgeno, tempfilels)
			msg = string("save imputed genofile in ",outfile,
				", tused=", round(tused,digits=1),"s")
			printconsole(io,verbose,msg)
		end		
	finally
		rm.(genofilels,force=true)		
		rm.(tempfilels,force=true)		
        rm(jltempdir,force=true,recursive=true)
	end
	tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, io, tused,"magicimpute_offspring!"; verbose)
	magicgeno
end

function save_imputed_gtgp(outfile::AbstractString,
	magicgeno::MagicGeno,tempfilels::AbstractVector;
	commentstring::AbstractString="##")
	outext = last(MagicBase.split_allext(basename(outfile)))
	outext in [".vcf",".vcf.gz"] || @error string("outfile ext must be in [.vcf,.vcf.gz]")
	isvcf = true
	delim = isvcf ? '\t' : ','
	myopen = outext in [".vcf.gz"] ? GZip.open : open
	myopen(outfile, "w") do io
		cc = commentstring
		write(io, cc*"fileformat=VCFv4.3\n")
		filedate = string(cc*"filedate=",string(now()))
		write(io, filedate,"\n")
		write(io, cc*"source=RABBIT/MagicImpute\n")
		msg = cc*"FORMAT=<ID=GP,Number=.,Type=Integer,Description=\"Probabilities of bi-allelic genotypes [0/0,1/1] (or [0/0,0/1,1/1] or [0/0,0/1,1/0,1/1]) for the vector length of 2 (or 3 or 4)\">"
		write(io, msg,"\n")
		msg = cc*"FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">"
		write(io, msg,"\n")
		msg = cc*"FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
		write(io, msg,"\n")
		msg = cc*"INFO=<ID=LINKAGEGROUP,Number=1,Type=String,Description=\"Linkage group in genetic map\">"
        write(io, msg,"\n")
		msg = cc*"INFO=<ID=POSCM,Number=1,Type=Float,Description=\"Genetic marker position in centiMorgan\">"
		write(io, msg,"\n")
		msg = cc*"INFO=<ID=FOUNDERERROR,Number=1,Type=Float,Description=\"Founder allelic error rate\">"
		write(io, msg,"\n")
		msg = cc*"INFO=<ID=OFFSPRINGERROR,Number=1,Type=Float,Description=\"Offspring allelic error rate\">"
		write(io, msg,"\n")
		msg = cc*"INFO=<ID=IMPUTETHRESH,Number=1,Type=Float,Description=\"Threshold of imputing offspring genotypes\">"
		write(io, msg,"\n")
		colnames = MagicBase.geno_colnames(magicgeno,isvcf,"all")
	    write(io, join(colnames,delim),"\n")
		nchr = length(tempfilels)
		# noff = size(magicgeno.magicped.offspringinfo,1)
		# progress = Progress(nchr;desc = "saving postprob..")
		for chr in 1:nchr			
			chrid = magicgeno.markermap[chr][1,:linkagegroup]
			mtx = jldopen(tempfilels[chr], "r") do file
				file["chr"] == chrid || @error string("inconsistent chr: ", file["chr"], ", ", chrid)
				file["postprob_mtx"]
		    end
			writedlm(io,mtx,delim)
	        flush(io)
			# next!(progress)
		end
	end
	outfile
end

function impute_offspring_chr!(genofile::AbstractString, 
    model::AbstractString,
    magicprior::NamedTuple;
	likeparam::LikeParam,    
	threshimpute::Real,
	israndallele::Bool,
    isfounderinbred::Bool,	
	phasealg::AbstractString,
	logio::Union{Nothing,AbstractString,IO},
    tempfile::AbstractString,
    verbose::Bool=true)
	io = isa(logio,AbstractString) ? open(logio,"w") : logio	
	try		
		startt = time()
		tsleep = parse(Int,replace(last(split(genofile,"_")),".csv.gz"=>"", "tsleep"=>""))
		sleep(tsleep)		
		tused = @elapsed magicgeno = readmagicgeno(genofile)
		chr = 1
	    chrid = magicgeno.markermap[chr][1,:linkagegroup]		
		mem1 = round(Int, memoryuse()/10^6)
		GC.gc()
		mem2 = round(Int, memoryuse()/10^6)
		printconsole(io,verbose,string("chr=", chrid, ", reading magicgeno, tused=", round(tused,digits=1),"s",
			", mem=",join([mem1,mem2],"|")))				
		# chrdecodefile saves HMM docode results for each offspring
		chrdecodefile =  string(tempname(dirname(tempfile); cleanup = false), "_decode_",chrid,".jld2")
		isphaseoffspring = !occursin(r"^unphase",phasealg)
		try	
			MagicReconstruct.reconstruct_chr!(magicgeno,chr, model,magicprior;
				likeparam,israndallele, isfounderinbred, 
				usepermarkererror = true, # using inferred error rate in magicgeno.markermap 
				hmmalg="forwardbackward", resetvirtual=true, 
				decodetempfile=chrdecodefile,logio=io,verbose
			)		
			if isfounderinbred
				chrfhaplo = magicgeno.foundergeno[chr]
			else
				chrfhaplo = permutedims(reduce(hcat,[reduce(vcat,i) for i in eachrow(magicgeno.foundergeno[chr])]))
			end			
			ischrx = magicgeno.markermap[chr][1,:linkagegroup] == "chrx"
			ischrx && error(string("TODO for sex chromosome"))
			# tempfile:  marker genoprob for each offspring			
			popmakeup = MagicReconstruct.calpopmakeup(magicgeno,chr, model,magicprior; isfounderinbred)			
			inputepsf = MagicBase.get_likeproperty(likeparam, :foundererror)			
			epsfls = MagicBase.get_errorls(magicgeno,chr;errorname = :foundererror, verbose=false)			
			if isnothing(epsfls)
				chrepsf = inputepsf
			else
				chrepsf = [ismissing(i) ? inputepsf : i for i in epsfls]				
			end
			condprob2postprob!(tempfile,chrdecodefile,chrfhaplo,chrepsf,popmakeup; isphaseoffspring, probdigits = 4)			
			if phasealg == "viterbi"				
				tviterbi = @elapsed MagicReconstruct.reconstruct_chr!(magicgeno,chr, model,magicprior;
					likeparam,israndallele, isfounderinbred, 
					usepermarkererror = true, hmmalg="viterbi", resetvirtual=true, 
					decodetempfile=chrdecodefile,logio=nothing,verbose = false
				)	
				printconsole(io,verbose,string("chr=", chrid, ", viterbi, tused=", round(tviterbi,digits=1),"s"))								
				viterbi2calledgeno!(magicgeno.offspringgeno[chr], chrdecodefile, chrfhaplo, chrepsf, popmakeup; threshimpute)				 
			else				
				callfrompostprob!(magicgeno.offspringgeno[chr],tempfile, threshimpute)				
			end
			offformat =  isphaseoffspring ? "GT_phased"  : "GT_unphased" 
			magicgeno.markermap[chr][!,:offspringformat] .= offformat				
			savemagicgeno(genofile,magicgeno)			
		finally			
			rm(chrdecodefile; force=true)
		end
		# tempfile: marker genoprob for each offspring -> prob_matrix (nsnp x noff)		
		gtgp2tempfile!(magicgeno,threshimpute,tempfile)				
		nsnp = size(magicgeno.markermap[chr],1)
		chrid = magicgeno.markermap[chr][1,:linkagegroup]
		mem1 = round(Int, memoryuse()/10^6)
		magicgeno.offspringgeno = nothing
		GC.gc()
		mem2 = round(Int,memoryuse()/10^6)
		msg = string("chr=", chrid, ", #marker=", nsnp, 
			", tused=", round(time()-startt,digits=1),"s",
			", mem=",mem1,"|",mem2, "MB, end",
		)
	    MagicBase.printconsole(io,verbose,msg)						
	finally
		isa(logio,AbstractString) && close(io)
	end
end

function gtgp2tempfile!(magicgeno::MagicGeno, threshimpute::Real,tempfile::AbstractString)	
	fformat = unique(magicgeno.markermap[1][!,:founderformat])
	fformat in [["GT_haplo"],["GT_phased"]] || @error string("unexpected founder format: ",fformat)
	calledmtx = MagicBase.togenomtx(magicgeno, 1; isvcf=true, target="all")	
	noff = size(magicgeno.offspringgeno[1],2)	
	jldopen(tempfile, "r") do file
	    for off in 1:noff
	        if haskey(file,string(off))
	            magicgeno.offspringgeno[1][:,off] .= file[string(off)]
	        else
	            @error string("postprob of offspring =",off," not exist")
	        end
	    end
	end	
	# togenomtx requires appropriate offspringformat!
	magicgeno.markermap[1][!,:offspringformat] .= "GP"	
	probmtx = MagicBase.togenomtx(magicgeno, 1; isvcf=true, target="all")			
	nf = size(magicgeno.magicped.founderinfo,1)
	gt = view(calledmtx,:,10+nf:size(calledmtx,2))
	gp = view(probmtx,:,10+nf:size(probmtx,2))	
	gp[1][1:2] == ".:" || @error string("unexpected offspring genotype ",gp[1])
	map!((x,y)->string(x,chop(y,head=1,tail=0)),gp,gt,gp)	
	# vcf colnames = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL","FILTER", "INFO", "FORMAT"]
	threshinfo = string("IMPUTETHRESH=",Float32(threshimpute))
	probmtx[:,8] .= [i=="." ? threshinfo : i*";"*threshinfo for i in probmtx[:,8]]
	chrid = magicgeno.markermap[1][1,:linkagegroup]
	jldopen(tempfile, "w") do file
		write(file,"chr",chrid)
		write(file, "postprob_mtx",probmtx)
	end
end

function callfrompostprob!(res::AbstractMatrix, tempfile::AbstractString,threshimpute::Real)	
	jldopen(tempfile,"r") do file		
		noffspring = file["noffspring"]
		nmarker = file["nmarker"]		
		size(res) == (nmarker, noffspring) || @error "inconsistent dimension"
		for off in 1:noffspring
			postprob = file[string(off)]
		    res[:,off] .= MagicBase.callfromprob.(postprob,threshimpute;isphased=true,ishaplo=false)			
		end		
	end
	res
end

function condprob2postprob!(outpostprobfile::AbstractString,		
	chrcondprobfile::AbstractString, 		
	chrfhaplo::AbstractMatrix,	
	chrepsf::Union{Real,AbstractVector}, 
	popmakeup::AbstractDict;	
	isphaseoffspring::Bool, 
	probdigits::Integer)		
	nstate, nfgl = MagicReconstruct.hmm_nstate_nfgl(popmakeup)    
	nsnp,nfgl2 = size(chrfhaplo)
	nfgl == nfgl2 || @error "inconsistent nfgl"		
	diplostates = MagicBase.prior_diploindex(nfgl) 	
	jldopen(chrcondprobfile,"r") do probfile				
		snporder = probfile["snporder"]
		length(snporder) == nsnp || @error "inconsistent #markers"			
		issetequal(snporder, 1:length(snporder)) || @error "inconsistent snps"			
		tindex = sortperm(snporder) # tindex[i]: t index for inputsnp index i			
		jldopen(outpostprobfile,"w") do file
			noffspring = sum([length(i["offspring"]) for i in values(popmakeup)])
			write(file, "noffspring", noffspring)			
			write(file, "nmarker", nsnp)						
			for popid in keys(popmakeup)
				offls = popmakeup[popid]["offspring"]
				nzstate =  popmakeup[popid]["nzstate"]
				nzorigin =  popmakeup[popid]["nzorigin"]				
				ishaploid =  popmakeup[popid]["ishaploid"]						
				if ishaploid
					if isa(chrepsf,AbstractVector)
						prior_haplo_epsf = [MagicReconstruct.haploprior(i) for i in chrepsf]
					else
						prior_haplo_epsf = MagicReconstruct.haploprior(chrepsf)
					end
					nzcol = MagicReconstruct.get_nzcol(nstate,popmakeup,popid)           		
					fderivehaplo = MagicReconstruct.calfderive(chrfhaplo; nzstate,ishaploid)								
					for off in offls						
						offprob = last(read(probfile, string(off)))'
						if isa(chrepsf,AbstractVector)
							postprob = [begin 
								priorsnp = prior_haplo_epsf[snp]									
								round.(priorsnp[:,fderivehaplo[snp,nzstate]] * offprob[nzcol,tindex[snp]],digits=probdigits) 
							end for snp in 1:nsnp]	
						else
							priorsnp = prior_haplo_epsf
							postprob = [round.(priorsnp[:,fderivehaplo[snp,nzstate]] * offprob[nzcol,tindex[snp]],digits=probdigits)  for snp in 1:nsnp]	
						end
						write(file, string(off), postprob)
					end
				else					
					nstate == nfgl^2 || error(string("mismatch between nfgl=",nfgl, " and nstate=",nstate))		
					if isa(chrepsf,AbstractVector)
						prior_diplo_epsf = [MagicReconstruct.diploprior(i) for i in chrepsf]
					else
						prior_diplo_epsf = MagicReconstruct.diploprior(chrepsf)
					end
					if isphaseoffspring							
						genostates = MagicBase.prior_genoindex(nfgl) 					
						b = [in(reverse(i),nzorigin) && !in(i,genostates) for i in nzorigin] 
						nzorigin_unphased = nzorigin[.!b]
						for i in eachindex(genostates)
							if !in(genostates[i], nzorigin_unphased) && in(reverse(genostates[i]), nzorigin_unphased) 
								reverse!(genostates[i])							
							end
						end
						issubset(nzorigin_unphased,genostates) || @error string("origins=",nzorigin_unphased, ", is not a subset of ancestral genotypes =",genostates)
						d2g = MagicBase.rulediplo2geno(diplostates, genostates; ismtx=true)
						ibdbool = allequal.(genostates)
						nonibdbool = .!ibdbool											
						b = [!in(i, nzorigin_unphased) for i in genostates]
						ibdbool[b] .= false
						nonibdbool[b] .= false	
						fderivegeno = calfderivegeno(chrfhaplo,genostates, nzorigin_unphased)						
						for off in offls							
							offprob = last(read(probfile, string(off))) * d2g  # transform diploporob into genoprob			
							offprob_ibd = view(offprob, :,ibdbool)'
							offprob_nonibd = view(offprob, :,nonibdbool)'
							if isa(chrepsf,AbstractVector)
								postprob = [begin								
									priorsnp = prior_diplo_epsf[snp]
									ls = priorsnp.ibd[:, fderivegeno[snp,ibdbool]] * offprob_ibd[:,tindex[snp]]
									ls .+= priorsnp.nonibd[:, fderivegeno[snp,nonibdbool]] * offprob_nonibd[:,tindex[snp]]							
									round.(ls,digits=probdigits)  
								end for snp in 1:nsnp]			
							else
								priorsnp = prior_diplo_epsf
								postprob = [begin																	
									ls = priorsnp.ibd[:, fderivegeno[snp,ibdbool]] * offprob_ibd[:,tindex[snp]]
									ls .+= priorsnp.nonibd[:, fderivegeno[snp,nonibdbool]] * offprob_nonibd[:,tindex[snp]]							
									round.(ls,digits=probdigits)  
								end for snp in 1:nsnp]			
							end
							write(file, string(off), postprob)
						end
					else							
						diplostates = MagicBase.prior_diploindex(nfgl) 		
						ibdbool = allequal.(diplostates)
						nonibdbool = .!ibdbool											
						b = [!in(i, nzorigin) for i in diplostates]
						ibdbool[b] .= false
						nonibdbool[b] .= false	
						fderivediplo = MagicReconstruct.calfderive(chrfhaplo; nzstate,ishaploid)										
						for off in offls							
							offprob = last(read(probfile, string(off)))
							offprob_ibd = view(offprob, :,ibdbool)'
							offprob_nonibd = view(offprob, :,nonibdbool)'
							if isa(chrepsf,AbstractVector)
								postprob = [begin
									priorsnp = prior_diplo_epsf[snp]	
									ls = priorsnp.ibd[:, fderivediplo[snp,ibdbool]] * offprob_ibd[:,tindex[snp]]
									ls .+= priorsnp.nonibd[:, fderivediplo[snp,nonibdbool]] * offprob_nonibd[:,tindex[snp]]							
									round.([ls[1],ls[2]+ls[3],ls[4]],digits=probdigits)  
								end for snp in 1:nsnp]			
							else
								priorsnp = prior_diplo_epsf
								postprob = [begin									
									ls = priorsnp.ibd[:, fderivediplo[snp,ibdbool]] * offprob_ibd[:,tindex[snp]]
									ls .+= priorsnp.nonibd[:, fderivediplo[snp,nonibdbool]] * offprob_nonibd[:,tindex[snp]]							
									round.([ls[1],ls[2]+ls[3],ls[4]],digits=probdigits)  
								end for snp in 1:nsnp]			
							end	
							write(file, string(off), postprob)
						end
					end
				end
			end			
		end		
	end
	outpostprobfile
end

function calfderivegeno(fhaplo::AbstractMatrix, genostates::AbstractVector, nzorigin_unphased::AbstractVector)    
    issubset(nzorigin_unphased,genostates) || @error string("unphased origins=",nzorigin_unphased, ", is not a subset of ancestral genotypes =",genostates)
    dict = Dict(genostates .=> 1:length(genostates))
    cols = [dict[i] for i in nzorigin_unphased]
    derrule=Dict(["NN", "N1", "1N", "N2", "2N", "11", "12", "21", "22"] .=> 1:9)
    nsnp=size(fhaplo,1)    	
	fderive=spzeros(Int,nsnp,length(genostates))        
    fderive[:,cols]=[get(derrule,join(fhaplo[snp,i]),-1) for snp=1:nsnp,i=nzorigin_unphased]
    if -1 in unique(fderive[:,cols])
        error("there exist unknown founder alleles")
    end
    fderive
end


function viterbi2calledgeno!(res::AbstractMatrix, 	
	chrdecodefile::AbstractString, 			
	chrfhaplo::AbstractMatrix, 
	chrepsf::Union{Real,AbstractVector}, 
	popmakeup::AbstractDict;		
	threshimpute::Real)
	nstate, nfgl = MagicReconstruct.hmm_nstate_nfgl(popmakeup)    
	nsnp,nfgl2 = size(chrfhaplo)
	nfgl == nfgl2 || @error "inconsistent nfgl"		
	size(res,1) == nsnp|| @error "inconsistent #markers"
	jldopen(chrdecodefile,"r") do viterbifile
		snporder = viterbifile["snporder"]
		length(snporder) == nsnp || @error "inconsistent #markers"
		issetequal(snporder, 1:length(snporder)) || @error "inconsistent snps"
		tindex = sortperm(snporder) # tindex[i]: t index for inputsnp index i				
		for popid in keys(popmakeup)
			offls = popmakeup[popid]["offspring"]
			nzstate =  popmakeup[popid]["nzstate"]				
			ishaploid =  popmakeup[popid]["ishaploid"]						
			if ishaploid
				if isa(chrepsf,AbstractVector)
					prior_haplo_epsf = [MagicReconstruct.haploprior(i) for i in chrepsf]
				else
					prior_haplo_epsf = MagicReconstruct.haploprior(chrepsf)
				end
				nzcol = MagicReconstruct.get_nzcol(nstate,popmakeup,popid)        
				iscolstate = nzcol == nzstate
				if !iscolstate
					nzcol2nzstate = Dict(nzcol .=> nzstate)
				end
				fderivehaplo = MagicReconstruct.calfderive(chrfhaplo; nzstate,ishaploid)							
				for off in offls						
					nzcol2, vpath = read(viterbifile,string(off))					
					nzcol2 == nzcol || @error "inconsistent nzcol"		
					if !iscolstate
						vpath .= [nzcol2nzstate[i] for i in vpath]
					end		
					if isa(chrepsf,AbstractVector)
						probls = [begin 
							priorsnp = prior_haplo_epsf[snp]
							beststate = vpath[tindex[snp]]
							derived = fderivehaplo[snp,beststate] 
							priorsnp[:,derived]
						end for snp in 1:nsnp]
					else
						priorsnp = prior_haplo_epsf
						probls = [begin 							
							beststate = vpath[tindex[snp]]
							derived = fderivehaplo[snp,beststate] 
							priorsnp[:,derived]
						end for snp in 1:nsnp]
					end
					res[:,off] .= MagicBase.callfromprob.(probls,threshimpute; ishaplo=false)			
				end
			else						
				nstate == nfgl^2 || error(string("mismatch between nfgl=",nfgl, " and nstate=",nstate))					
				if isa(chrepsf,AbstractVector)
					prior_diplo_epsf = [MagicReconstruct.diploprior(i) for i in chrepsf]
				else
					prior_diplo_epsf = MagicReconstruct.diploprior(chrepsf)
				end
				diplostates = MagicBase.prior_diploindex(nfgl)
				ibdbool = allequal.(diplostates)					
				fderivediplo = MagicReconstruct.calfderive(chrfhaplo; nzstate,ishaploid)						
				for off in offls					
					nzcol, vpath = read(viterbifile,string(off))					
					nzcol == nzstate || @error "inconsistent nzcol and nzstate"
					if isa(chrepsf,AbstractVector)
						probls = [begin 
							priorsnp = prior_diplo_epsf[snp]
							beststate = vpath[tindex[snp]]
							derived = fderivediplo[snp,beststate] 
							isibd = ibdbool[beststate]
							isibd ? priorsnp.ibd[:,derived] : priorsnp.nonibd[:,derived] 
						end for snp in 1:nsnp]
					else
						priorsnp = prior_diplo_epsf
						probls = [begin 							
							beststate = vpath[tindex[snp]]
							derived = fderivediplo[snp,beststate] 
							isibd = ibdbool[beststate]
							isibd ? priorsnp.ibd[:,derived] : priorsnp.nonibd[:,derived] 
						end for snp in 1:nsnp]
					end
					res[:,off] .= MagicBase.callfromprob.(probls,threshimpute;isphased=true,ishaplo=false)			
				end
			end
		end				
	end
	res
end

