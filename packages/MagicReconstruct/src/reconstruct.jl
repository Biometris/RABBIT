
function reconstruct!(magicgeno::MagicGeno;
    model::AbstractString="jointmodel",
	likeparameters::LikeParameters=LikeParameters(),   
	israndallele::Bool=true, 
	isfounderinbred::Bool=true,		
    hmmalg::AbstractString="forwardbackward",
	posteriordigits::Integer = 4, 
    isMMA::Bool=false,
    isparallel::Bool=false,		
    outstem::Union{Nothing,AbstractString}="outstem",    	
	outext::AbstractString=".csv.gz",
	workdir::AbstractString=pwd(),
    tempdirectory::AbstractString = tempdir(),	
    io::Union{Nothing,IO}=nothing,
    verbose::Bool=true)		
	magicprior=MagicReconstruct.calmagicprior(magicgeno,model; isfounderinbred)	
	starttime = time()
    nchr = length(magicgeno.markermap)
    if isMMA
        reconstructfun = reconstruct_chr_mma
    else
        reconstructfun = reconstruct_chr!
    end
	nsnpls = [size(i,1) for i in magicgeno.markermap]
	chroo = MagicBase.get_chroo(nsnpls)		
    jltempdir = mktempdir(tempdirectory; prefix="jl_magicreconstruct_", cleanup=true)
	tempid = tempname(jltempdir,cleanup=false)
    decodefilels = [string(tempid, "_reconstruct_chr",chr,".tmp") for chr in 1:nchr]	
	# split and save geno by chromosome
	tempid = tempname(jltempdir,cleanup=false)
	tused = @elapsed genofilels = MagicBase.saveby_chromosome(magicgeno; nworker = nworkers(),
		workdir=jltempdir, outstem=tempid*"_submagicgeno")
	mem1 = round(Int, memoryuse()/10^6)
	magicgeno.markermap = nothing
	magicgeno.foundergeno = nothing # foundergeno is removed to save memory
	magicgeno.offspringgeno = nothing # offgenogeno is removed to save memory
	GC.gc()
	GC.gc()
	mem2 = round(Int, memoryuse()/10^6)
	msg = string("saveby_chromosome, tused=",round(tused,digits=1), 
		"s, mem=",mem1,"|",mem2,"MB")
	printconsole(io,verbose,msg)	
	# keep consistent with MagicBase.savemagicancestry_targz, so that it can be read by MagicBase.readmagicancestry_targz
	resdir_ancestry = mktempdir(workdir; prefix="jl_savemagicancestry_", cleanup=true)	
	resfilels = [begin 
		chrlab = MagicBase.get_ancestry_chrlab(chr,nchr)
    	abspath(resdir_ancestry,string(outstem,"_",chrlab, outext))
	end for chr in 1:nchr]	
    try
        if isparallel && nprocs()>1			
			tempid = tempname(jltempdir,cleanup=false)
			logfilels = [string(tempid, "_chr",chr,".log") for chr in 1:nchr]	            			
			try
				pmap((x,y,z,w)->reconstructfun(x, model,magicprior;
					likeparameters, israndallele, isfounderinbred, hmmalg,posteriordigits,				
					decodetempfile=y, chroutfile=z, logio=w,verbose),
					genofilels[chroo],decodefilels[chroo],resfilels[chroo],logfilels[chroo])
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
			end
			rm.(logfilels,force=true)
			@everywhere GC.gc()
        else
            for chr in eachindex(genofilels)
				reconstructfun(genofilels[chr], model,magicprior;
					likeparameters, israndallele, isfounderinbred, hmmalg, posteriordigits,				
					decodetempfile=decodefilels[chr], chroutfile=resfilels[chr], logio=io,verbose) 
			end
        end        		        
		msg = string("tused =", round(time()-starttime,digits=1), "s for all chromosomes")
		MagicBase.printconsole(io,verbose,msg)				        
		if outext == ".jld2"			
			outfile = getabsfile(workdir, string(outstem,"_ancestry.jld2.tar.gz"))
    		MagicBase.create_targz(resdir_ancestry,outfile)
		else			
			outfile = getabsfile(workdir, string(outstem,"_ancestry.csv.gz.tar"))
    		MagicBase.create_tar(resdir_ancestry,outfile)
		end		
		outfile
    finally
		rm.(genofilels,force=true)		
        rm.(decodefilels,force=true)		
        rm(jltempdir,force=true,recursive=true)
		rm.(resfilels,force=true)		
		rm(resdir_ancestry,force=true,recursive=true)		
    end
end

function reconstruct_chr!(genofile,
    model::AbstractString,	
    magicprior::NamedTuple;
	likeparameters::LikeParameters=LikeParameters(),   
	israndallele::Bool=true, 
	isfounderinbred::Bool=true,	
    hmmalg::AbstractString="forwardbackward",
	posteriordigits::Integer = 4, 
    logio::Union{Nothing,AbstractString,IO},
    decodetempfile::AbstractString,
	chroutfile::AbstractString, 
    verbose::Bool=true)
	tsleep = parse(Int,replace(last(split(genofile,"_")),".csv.gz"=>"", "tsleep"=>""))
	sleep(tsleep)	
	tused = @elapsed magicgeno = readmagicgeno(genofile)
	chr = 1
	chrid = magicgeno.markermap[chr][1,:linkagegroup]		
	mem1 = round(Int, memoryuse()/10^6)
	GC.gc()
	mem2 = round(Int, memoryuse()/10^6)
	io = isa(logio, AbstractString) ? nothing : logio
	printconsole(io,verbose,string("chr=", chrid, ", reading magicgeno, tused=", round(tused,digits=1),"s",
		", mem=",join([mem1,mem2],"|")))			
	chrres = reconstruct_chr!(magicgeno, chr,  model,magicprior;
		likeparameters,israndallele, isfounderinbred,hmmalg,posteriordigits,
		decodetempfile,logio,verbose) 		
	chrancestry = combine2ancestry(magicgeno,chr, model, chrres, isfounderinbred,posteriordigits)
	outext = last(MagicBase.split_allext(chroutfile))
	if outext == ".jld2"				
		JLD2.save(chroutfile, chrid, chrancestry)
	elseif outext == ".csv.gz"
		savemagicancestry(chroutfile,chrancestry)		
	else
		error(string("unknown outext=",outext))
	end
	chroutfile
end

function reconstruct_chr_mma(genofile,
    model::AbstractString,
    magicprior::NamedTuple;
	likeparameters::LikeParameters=LikeParameters(),   
	israndallele::Bool, 
	isfounderinbred::Bool=true,	
    hmmalg::AbstractString="forwardbackward",
	posteriordigits::Integer = 4, 
    logio::Union{Nothing,AbstractString,IO},
    decodetempfile::AbstractString,
	chroutfile::AbstractString, 
    verbose::Bool=true)
	tsleep = parse(Int,replace(last(split(genofile,"_")),".csv.gz"=>"", "tsleep"=>""))
	sleep(tsleep)	
	tused = @elapsed magicgeno = readmagicgeno(genofile)
	chr = 1
	chrid = magicgeno.markermap[chr][1,:linkagegroup]		
	mem1 = round(Int, memoryuse()/10^6)
	GC.gc()
	mem2 = round(Int, memoryuse()/10^6)
	io = isa(logio, AbstractString) ? nothing : logio
	printconsole(io,verbose,string("chr=", chrid, ", reading magicgeno, tused=", round(tused,digits=1),"s",
		", mem=",join([mem1,mem2],"|")))				
		chrres = reconstruct_chr_mma(magicgeno, chr,  model,magicprior;
		likeparameters,israndallele, isfounderinbred,hmmalg,posteriordigits,
		decodetempfile,logio,verbose) 		
	chrancestry = combine2ancestry(magicgeno,chr, model, chrres, isfounderinbred,posteriordigits)
	savemagicancestry(chroutfile,chrancestry)
	chroutfile
end


function combine2ancestry(magicgeno::MagicGeno, chr::Integer, model::AbstractString, 
	chrres::Tuple,isfounderinbred::Bool,posteriordigits::Integer)
    chrloglike,chrdecodefile = chrres
	hmmalg = load(chrdecodefile,"hmmalg")
    loglike = reshape(chrloglike,1, :)
    decode = [[j[2] for j in last(MagicReconstruct.read_decodefile(chrdecodefile))]]    
    if hmmalg == "forwardbackward"
        if model=="depmodel"
            haploprob=decode
            genoprob=diploprob=nothing
        else
            diploprob=decode
            genoprob=haploprob=nothing
        end
        viterbipath = nothing
    else
        # hmmalg == "viterbi"
        viterbipath = [reduce(hcat,i) for i in decode]
        diploprob = genoprob = haploprob = nothing        
    end
    statespace = MagicReconstruct.getpriorstatespace(magicgeno.magicped; isfounderinbred)
    magicancestry=MagicAncestry(magicgeno.magicped,magicgeno.markermap[[chr]],magicgeno.foundergeno[[chr]],
        statespace, viterbipath,diploprob,genoprob,haploprob,nothing,loglike,magicgeno.misc)
    if model != "depmodel" 
		MagicBase.setinbredcoef!(magicancestry; posteriordigits)
		if hmmalg == "forwardbackward"
			MagicBase.setgenoprob!(magicancestry; posteriordigits)
			MagicBase.sethaploprob!(magicancestry; posteriordigits)
		end
    end
    magicancestry
end

function reconstruct_chr_mma(magicgeno::MagicGeno,chr::Integer,
    model::AbstractString,
    magicprior::NamedTuple;
	likeparameters::LikeParameters=LikeParameters(),   	
	israndallele::Bool=true, 
    isfounderinbred::Bool=true,	
    hmmalg::AbstractString="forwardbackward",
	posteriordigits::Integer = 4, 
    logio::Union{Nothing,AbstractString,IO},
    decodetempfile::AbstractString,
    verbose::Bool=true)
	io = isa(logio,AbstractString) ? open(logio,"w") : logio
	try
		starttime = time()
		_, epsf, epso, epso_perind, seqerror, allelebalancemean,allelebalancedisperse,alleledropout = MagicBase.extract_likeparameters(likeparameters)		
	    # inter-marker distance in Morgan
	    deltd= diff(magicgeno.markermap[chr][!,:poscm]) * 0.01
		isoffphased = false
	    isautosome, issnpGT, ishaploidls, fderive, offcode =
	        precompute_chr_mma(magicgeno,chr,model,isoffphased, isfounderinbred)
	    prior = isautosome ? magicprior.autosome : magicprior.allosome
	    nfgl=size(magicgeno.foundergeno[chr],2)*(isfounderinbred ? 1 : 2)
	    nstate= (false in ishaploidls) ? nfgl^2 : nfgl
	    nsnp,noff=size(magicgeno.offspringgeno[chr])
	    loglike=Vector{Float64}(undef,noff)
	    subpopls=collect(keys(prior))		
	    jldopen(decodetempfile,"w") do file			
	        for subpopid in subpopls
	            subpop=findall(magicgeno.magicped.offspringinfo[!,:member] .== subpopid)
	            ishaploidsubls=ishaploidls[subpop]
	            if length(union(ishaploidsubls)) !=1
	                @error string("both diploid and haploid are include in subpopulation ",subpopid)
	            end
	            # calculate initprob and tranprobseq; take only accessible states
	            # nzindex refers to the set of indices with states being acessible
	            nzindex,initprob,tranrate = prior[subpopid]
	            # abs is to assure positive tranprob; negative due to numerical error
	            tranprobseq=[get_tranprobmatrix(i,tranrate) for i = deltd]
	            ishaploidsub = all(ishaploidsubls)
	            if ishaploidsub && (false in ishaploidls)
	                # to re-write for male x chromsome
	                dict=Dict([(i-1)*nfgl+i=>i for i=1:nfgl])
	                datanzindex=[get(dict,i,0) for i=nzindex]
	            else
	                datanzindex=nzindex
	            end				
				dataprobseq = [zeros(_float_like,length(nzindex))  for _ in 1:nsnp]
	            for off in subpop
	                obsseq = offcode[:,off]					
					caldataprobseq!(dataprobseq,obsseq,epsf,epso,epso_perind, seqerror, allelebalancemean,allelebalancedisperse,alleledropout,
                    	fderive,nzstate,isoffphased,israndallele,issnpGT,ishaploidsub)
	                dataprobseq2 = [Vector(dataprobseq[snp][datanzindex]) for snp=1:nsnp]					
	                if hmmalg == "forwardbackward"
	                    loglike[off],prob= HMM.posteriordecode(initprob,tranprobseq,dataprobseq2;digits=posteriordigits)
	                    off_decode = spzeros(nsnp,nstate)
	                    off_decode[:,nzindex] .= round.(reduce(hcat,prob)',digits=4)
	                else
	                    # hmmalg == "viterbi"
	                    off_decode = spzeros(nsnp)
	                    loglike[off],vpath = HMM.viterbidecode(initprob,tranprobseq,dataprobseq2)
	                    off_decode[nzindex] .= nzindex[vpath]
	                end
	                write(file, string(off), off_decode)
	            end
	        end
	    end
	    chrid = magicgeno.markermap[chr][1,:linkagegroup]
	    nsnp = size(magicgeno.markermap[chr],1)
	    msg = string("chr=", chrid, ", #snp=", nsnp, ", tused=",
	        round(time()-starttime,digits=1)," in reconstruct")
	    MagicBase.printconsole(io,verbose,msg)
		(loglike, decodetempfile)
	finally
		isa(logio,AbstractString) && close(io)
	end	
end

function get_all_errls(magicgeno::MagicGeno,chr::Integer; 
	likeparameters::LikeParameters=LikeParameters(),   	
	io::Union{Nothing,IO}=nothing, verbose::Bool=true)
	# set error rates	
	isanyGP = in("GP", magicgeno.markermap[chr][!,:offspringformat])
	if isanyGP
		epsfls0 = nothing
		epsols0 = nothing
		epsols0_perind = nothing
	else	
		epsfls0 = MagicBase.get_errorls(magicgeno,chr;errorname = :foundererror, io,verbose)
		epsols0 = MagicBase.get_errorls(magicgeno,chr;errorname = :offspringerror, io,verbose)			
		epsols0_perind = MagicBase.get_errorls(magicgeno,chr;errorname = :peroffspringerror, io,verbose)		
	end
	isanyAD = in("AD", magicgeno.markermap[chr][!,:offspringformat])	
	if isanyAD
		seqerrorls0 = MagicBase.get_errorls(magicgeno,chr;errorname = :seqerror, io,verbose)
		allelebalancemeanls0 = MagicBase.get_errorls(magicgeno,chr;errorname = :allelebalancemean, io,verbose)
		allelebalancedispersels0 = MagicBase.get_errorls(magicgeno,chr;errorname = :allelebalancedisperse, io,verbose)
		alleledropoutls0 = MagicBase.get_errorls(magicgeno,chr;errorname = :alleledropout, io,verbose)	
	else
		seqerrorls0 = nothing
		allelebalancemeanls0 = nothing
		allelebalancedispersels0 = nothing
		alleledropoutls0 = nothing	
	end
	liketargetls, epsf, epso, epso_perind, seqerror, allelebalancemean,allelebalancedisperse,alleledropout = MagicBase.extract_likeparameters(likeparameters)		
	nsnp  = size(magicgeno.markermap[chr],1)	
	epsfls = merge_errorls(epsfls0,epsf,nsnp)	
	epsols = merge_errorls(epsols0,epso,nsnp)	
	noff = size(magicgeno.magicped.offspringinfo,1)
	epsols_perind = merge_errorls(epsols0_perind,epso_perind,noff)	
	seqerrorls = merge_errorls(seqerrorls0,seqerror,nsnp)	
	allelebalancemeanls = merge_errorls(allelebalancemeanls0,allelebalancemean,nsnp)		
	allelebalancedispersels = merge_errorls(allelebalancedispersels0,allelebalancedisperse,nsnp)		
	alleledropoutls = merge_errorls(alleledropoutls0,alleledropout,nsnp)	
	(liketargetls, epsfls, epsols, epsols_perind, seqerrorls,allelebalancemeanls,allelebalancedispersels,alleledropoutls)
end

function reconstruct_chr!(magicgeno::MagicGeno,chr::Integer,
    model::AbstractString,
    magicprior::NamedTuple;
	likeparameters::LikeParameters=LikeParameters(),   	
	usepermarkererror::Bool=false, 
	israndallele::Bool=true, 
	isfounderinbred::Bool=true,	
    hmmalg::AbstractString="forwardbackward",
	posteriordigits::Integer = 4, 
    logio::Union{Nothing,AbstractString,IO}=nothing,
    decodetempfile::AbstractString,
    verbose::Bool=true)
	io = isa(logio,AbstractString) ? open(logio,"w") : logio
	try		
		chrid = magicgeno.markermap[chr][1,:linkagegroup]				
		if isfounderinbred
			chrfhaplo = magicgeno.foundergeno[chr]
		else
			chrfhaplo = permutedims(reduce(hcat,[reduce(vcat,i) for i in eachrow(magicgeno.foundergeno[chr])]))
		end
		nsnp = size(chrfhaplo,1)
	    chroffgeno = magicgeno.offspringgeno[chr]
		issnpGT = [occursin("GT",i) for i in magicgeno.markermap[chr][!,:offspringformat]]		
		gtls = unique(magicgeno.markermap[chr][issnpGT,:offspringformat])
		isoffphased = gtls== ["GT_phased"]		
		length(gtls) > 1 && @error string("mixed phased and unphased genotypes, offspring format=",gtls)		
		popmakeup,priorprocess= calpriorprocess(magicgeno,chr, model,magicprior; isfounderinbred)				
		if usepermarkererror		
			_, epsf, epso, epso_perind, seqerror, allelebalancemean,allelebalancedisperse,alleledropout = get_all_errls(magicgeno,chr; likeparameters,io,verbose)					
		else
			_, epsf, epso, epso_perind, seqerror, allelebalancemean,allelebalancedisperse,alleledropout = MagicBase.extract_likeparameters(likeparameters)					
		end
		# decoding		
	    tused = @elapsed res = hmmdecode_chr(chrfhaplo,chroffgeno,popmakeup,priorprocess;
	        epsf,epso, epso_perind,seqerror, 
			allelebalancemean, allelebalancedisperse, alleledropout,  
			hmmalg, posteriordigits, issnpGT, isoffphased, israndallele, decodetempfile)				
		chrlen = round(Int,magicgeno.markermap[chr][end,:poscm] - magicgeno.markermap[chr][1,:poscm])
		mem1 = round(Int, memoryuse()/10^6)
		GC.gc()
		mem2 = round(Int,memoryuse()/10^6)
		msg = string("chr=", chrid, ", ", hmmalg, 				
			", #snp=", nsnp,
			", chrlen=", chrlen, "cM",
			", tused=",round(tused,digits=1),"s",		
			", mem=",mem1,"|",mem2, "MB",		
		)
		MagicBase.printconsole(io,verbose,msg)	
        res
	finally
		isa(logio,AbstractString) && close(io)
	end
end

function get_outerlier_indices(ls::AbstractVector; 
	tukeyfence::Real=3,	
	minoutlier=0.0,
	side::String="both")		
    q1,q3=quantile(ls,[0.25,0.75])
	if side in ["lower","both"]
		lowbound= q3-tukeyfence*(q3-q1)
		res = findall(ls .< lowbound)
	else
		res = Int[]
	end
	if side in ["upper","both"]
		upbound= max(minoutlier,q3+tukeyfence*(q3-q1))
		append!(res,findall(ls .> upbound))
	end
	res
end

function del_error_indices(errorls::AbstractVector;
	minoutlier=0.05,tukeyfence::Real=3,max_error::Real=0.2)	
	ls = log.(errorls)
    q1,q3=quantile(ls,[0.25,0.75])
	upbound=q3+tukeyfence*(q3-q1)	    
	upbound=min(max(upbound,log(minoutlier)),log(max_error))
    findall(ls .> upbound)
end

function merge_errorls(errls::Union{Nothing,AbstractVector},
	inputerr::Union{Real,AbstractVector},
	nsnp::Integer)
	if isnothing(errls)
		res = isa(inputerr,AbstractVector) ? inputerr : inputerr*ones(nsnp) 
	else
		nsnp == length(errls) || @error "inconsistent #markers"
		if isa(inputerr,AbstractVector)
			# @warn string("use epsf from markermap rather than input argment")
			res = map((x,y) -> ismissing(x) ? y : x, errls, inputerr)
		else
			res = [ismissing(i) ? inputerr : i for i in errls]
		end
	end
	res
end
