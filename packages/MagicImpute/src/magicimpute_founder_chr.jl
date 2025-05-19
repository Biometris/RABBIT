function impute_refine_repeat_chr!(magicgenofile::AbstractString,nrepeatimpute::Tuple; 
	magicprior::NamedTuple,
    model::AbstractString="jointmodel",
	israndallele::Bool,
	isfounderinbred::Bool=true,	        	
	byfounder::Integer=0, 
	startbyhalf::Integer,
	iscorrectfounder::Bool = true,	
	isimputefounder::Union{Nothing,Bool}=nothing, 
	isallowmissing::Bool, 
	threshproposal::Real,
    isdelmarker::Bool = true,
    delsiglevel::Real = 0.01,
	isinferjunc::Bool = false, 
	isinfererror::Bool = false,	    
	likeparameters::LikeParameters, 
	threshlikeparameters::ThreshLikeParameters, 
	priorlikeparameters::PriorLikeParameters,
	tukeyfence::Real=3.0,			
	minoutlier::Real=0.05, 
	isspacemarker::Bool = false,
    trimcm::Real=20,
	trimfraction::Real=0.05,  #cM
	skeletonsize::Union{Nothing,Integer},
	isordermarker::Bool = false,
	inputneighbor::Union{Nothing,AbstractDict}=nothing,
	slidewin::Union{Nothing,Integer} = nothing,  
	slidewin_neighbor::Union{Nothing,Integer} = 200,      
	orderactions::AbstractVector = ["inverse","permute"],  
    orderactions_neighbor::AbstractVector = ["inverse11","inverse01"],  
	inittemperature::Real= isordermarker ? 2.0 : 0.0,
    coolrate::Real=0.85,
    minaccept::Real=0.15,
	spacebyviterbi::Bool=false, 
    logio::Union{Nothing,AbstractString,IO},    		
	maxiter::Integer = 50,		
	imputetempfile::AbstractString, 
    verbose::Bool,	
	more_verbose::Bool)	
	io = isa(logio,AbstractString) ? open(logio,"w") : logio	
	nrepeatmin, nrepeatmax = nrepeatimpute
	try	
		chrid = last(split(first(splitext(basename(imputetempfile))),"_"))		
		tempdirectory = dirname(imputetempfile)
		tempid = tempname(tempdirectory,cleanup=true)		
		runtempfile = string(tempid, "_impute_founder_adptrepeat_",chrid, ".temp") 
		magicgenofilels = get_magicgenofilels(magicgenofile,nrepeatmax)
		for i in eachindex(magicgenofilels)
			cp(magicgenofile,magicgenofilels[i]; force=true)
		end		
		# get connected components of the populaiton
		magicped = readmagicped(magicgenofile) # magicgenofile contains both genodata and pedinfor
		pedcompoents = get_maigcped_connetedcomponents(magicped; isfounderinbred, model, magicprior)
		ncomponents = length(pedcompoents)
		# initialization
		genodiff = [-1*ones(nrepeatmax,nrepeatmax) for _ in 1:ncomponents]
		phasediff = [-1*ones(nrepeatmax,nrepeatmax) for _ in 1:ncomponents]
		tempio = IOBuffer()		
		try 		
			jldopen(runtempfile,"w+") do file			
				loglhis = Vector{Vector{Float64}}()
				bestrun_components = ones(Int, ncomponents)
				minrun_components = ones(Int, ncomponents)
				bestlogl_components = -Inf*ones(ncomponents)		
				isdone_components = falses(ncomponents)			
				genodiff_components = Any[NaN for _ in 1:ncomponents]				
				phasediff_components = Any[NaN for _ in 1:ncomponents]												
				for runit=1:nrepeatmax			
					# magicgenofilels[runit] saves calculated results
					# if isfounderinbred, founder genotypes in fhaplodf are alleles, and otherwise phased genotypes
					loglikels, fmissls, fhaplodf = impute_refine_chr!(magicgenofilels[runit]; 
						magicprior,
						model,israndallele, isfounderinbred,
						byfounder,startbyhalf, 
						isinferjunc, iscorrectfounder,isimputefounder,isallowmissing, threshproposal,
						isdelmarker, delsiglevel,
						isspacemarker, trimcm, trimfraction,skeletonsize,
						isinfererror, likeparameters, threshlikeparameters, priorlikeparameters, tukeyfence, minoutlier, 															
						isordermarker, inputneighbor,slidewin,slidewin_neighbor,orderactions,orderactions_neighbor,
						inittemperature, coolrate, minaccept, spacebyviterbi,                   
						logio=io, imputetempfile, repeatrun=runit, maxiter, verbose,more_verbose)																	
					write(file, string("run",runit),fhaplodf)
					logl_components = [sum(loglikels[pedcomponent["offspring"]]) for pedcomponent in pedcompoents]
					logl_components .= round.(logl_components,digits=1)													
					push!(loglhis, logl_components)		
					done_components = []  # collect compoents that converge after current run
					for c in findall(.!isdone_components)
						# update logl
						if logl_components[c] >= bestlogl_components[c]
							bestrun_components[c] = runit
							bestlogl_components[c] = logl_components[c]
						end
						bestrun = bestrun_components[c]
						founders = pedcompoents[c]["founders"]											
						bestfhaplo = file[string("run",bestrun)]					
						# compare runs by fhaplo
						if !isdone_components[c] 						
							genodiff[c][runit,runit]  = 0.0					
							phasediff[c][runit,runit]  = 0.0	
							if runit > 1								
								if bestrun == runit 										
									for run in 1:runit-1										
										runfhaplo = file[string("run",run)]
										genodiff[c][bestrun,run] = genodiff[c][run,bestrun] = caldiff_fhaplo(bestfhaplo,runfhaplo,founders; isfounderinbred)							 							
									end
									genodiff[c][runit,runit]  = 0.0
								else									
									genodiff[c][bestrun,runit] = genodiff[c][runit,bestrun] = caldiff_fhaplo(bestfhaplo,fhaplodf,founders; isfounderinbred)							 
								end
								if !isfounderinbred
									if bestrun == runit 										
										for run in 1:runit-1										
											runfhaplo = file[string("run",run)]
											phasediff[c][bestrun,run] = phasediff[c][run,bestrun] = calphasediff_fhaplo(bestfhaplo,runfhaplo,founders)							 							
										end
										phasediff[c][runit,runit]  = 0.0
									else									
										phasediff[c][bestrun,runit] = phasediff[c][runit,bestrun] = calphasediff_fhaplo(bestfhaplo,fhaplodf,founders)							 
									end									
								end
							end
							# comparse fhaplo between runs
							if runit >= nrepeatmin	|| runit == nrepeatmax				
								genodiff_best = round.(genodiff[c][bestrun,1:runit],digits=3)
								genodiff_components[c] = copy(genodiff_best)
								genodiff_best[bestrun] = Inf				
								minrunls = setdiff(findall(isinf.(genodiff_best)),[bestrun])		
								if isempty(minrunls)
									mindiff = minimum(genodiff_best)		
									minrun = rand(findall(genodiff_best .== mindiff))
								else
									minrun = first(minrunls)
									mindiff = 0.0
								end																
								minrun_components[c] = minrun
								# diff betetween best run and other runs
								genodiff_components[c] = round.(genodiff[c][bestrun,1:runit],digits=3)
								if !isfounderinbred																		
									phasediff_components[c] = round.(phasediff[c][bestrun,1:runit],digits=3)																		
								end											
								isapproxrun = mindiff <= 0.001 || abs(loglhis[minrun][c] - bestlogl_components[c]) <= log(10.0)
								if runit == nrepeatmax && !isapproxrun
									warnmsg = string("chr=",chrid,", run=", runit, ", optimization reaches nrepeatmax=",nrepeatmax, " for component=",c)
									verbose && @warn warnmsg
									printconsole(io, false, "Warning: "*warnmsg)
									write(tempio, warnmsg, "\n")
								end																									
								if runit == nrepeatmax || isapproxrun
									isdone_components[c] = true																						
									push!(done_components,c)											
									if mindiff > 0 
										minfhaplo = file[string("run",minrun)]
										isdiff = isdiff_fhaplo(bestfhaplo,minfhaplo,founders; isfounderinbred)		
										printconsole(io, verbose, string("chr=", chrid, ", #diff=",sum(isdiff), " between bestrun ", bestrun, " and neighborrun ", minrun))										
									end									
								end						
							end  
						end								
					end			
					if ncomponents == 1 && length(done_components) == 1
						c = only(done_components) 
						(all(isdone_components) && (c == 1)) || @error string("unexpected done_components=",done_components, ", or unexpected isdone_components=",isdone_components)						
						cp(magicgenofilels[bestrun_components[c]],magicgenofile; force=true)
					else
						if !isempty(done_components)
							# copy results of bestrun to magicgenofile and future-input genofiles
							magicgeno = readmagicgeno(magicgenofile)
							inputmarkermap = only(magicgeno.markermap)						
							markerdict = Dict(inputmarkermap[!,:marker] .=> 1:size(inputmarkermap,1))						
							for c in done_components
								founders = pedcompoents[c]["founders"]											
								bestgeno = readmagicgeno(magicgenofilels[bestrun_components[c]])
								posls = [markerdict[i] for i in only(bestgeno.markermap)[!,:marker]]						
								magicgeno.foundergeno[1][posls, founders] .= bestgeno.foundergeno[1][:, founders] 							
								inputmarkermap[posls,:founderformat] .= isfounderinbred ? "GT_haplo" : "GT"
							end										
							savemagicgeno(magicgenofile,magicgeno)																	
							# all the subsequenct runs start from the new estimated founder genotypes
							for i in runit+1:nrepeatmax				
								cp(magicgenofile,magicgenofilels[i]; force=true)
							end							
						end						
					end
					# print and check if done
					msg = string("chr=",chrid,", run=", runit)
					if ncomponents == 1
						msg *= string(", logl=",only(logl_components), 
							", bestrun=", only(bestrun_components), 							
							", bestlogl=",only(bestlogl_components),
							", loglhis=",only.(loglhis))
						if isa(only(genodiff_components),AbstractVector)
							msg *= string(", genodiff=",only(genodiff_components))
						end
						if !isfounderinbred
							if isa(only(phasediff_components),AbstractVector)
								msg *= string(", phasediff=",only(phasediff_components))
							end							
						end				
						if runit >= nrepeatmin > 1	
							msg *= string(", neighborrun=", only(minrun_components))								
						end
					else
						msg *= string(", logl=",logl_components, 
							", bestrun=", bestrun_components, 							
							", bestlogl=",bestlogl_components,							
							", loglhis=",loglhis)												
						if any(isa.(genodiff_components,AbstractVector))
							msg *= string(", genodiff=",genodiff_components)
						end
						if !isfounderinbred
							if any(isa.(phasediff_components,AbstractVector))
								msg *= string(", phasediff=",phasediff_components)
							end
						end						
						if runit >= nrepeatmin > 1
							msg *= string(", neighborrun=", minrun_components)							
						end
					end
					printconsole(io,verbose,msg)
					write(tempio, msg,"\n")
					if runit >= nrepeatmin && all(isdone_components)																		
						break
					end					
				end			
				all(isdone_components) || @error string("unfinished compoents = ",findall(.!isdone_components))
				isfounderinbred || resolve_GT!(magicgenofile; io,verbose)				
				msg = String(take!(tempio))
				msg = string("chr=",chrid, ", run_summary:\n", rstrip(msg,'\n'))
				printconsole(io,verbose,msg)
				nothing			
			end
		finally
			close(tempio)			
			rm(runtempfile; force=true)
			rm.(magicgenofilels; force=true) 
		end
	finally
		isa(logio,AbstractString) && close(io)		
	end
end

function resolve_GT!(magicgenofile::AbstractString; io::Union{IO, Nothing}=nothing, verbose=true) 
	# only for outbred parents
	# delete markers with "GT_unphased" for outbred parents, or to consider mixed phased/unphased markers
	magicgeno = readmagicgeno(magicgenofile)	
	for chr in eachindex(magicgeno.markermap)
		chrmarkermap = magicgeno.markermap[chr]				
		b = [in(i, ["GT_phased","GT"]) for i in chrmarkermap[!,:founderformat]]
		all(b) && continue		
		keepat!(chrmarkermap,b)
		magicgeno.foundergeno[chr] = magicgeno.foundergeno[chr][b,:]
		magicgeno.offspringgeno[chr] = magicgeno.offspringgeno[chr][b,:]
	end
	# at a given marker, some parents' genotypes may be unphased! set them to phased genotypes if they are homozygoutes and missing if heterozygous			
	nhalf = 0
	nhetero = 0
	for chr in eachindex(magicgeno.markermap)
		chrmarkermap = magicgeno.markermap[chr]				
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
						if in("N", gphased) 
							nhalf += 1
						else
							nhetero += 1
						end
					end								
				end
			end			
		end		
		chrmarkermap[!,:founderformat] .= "GT_phased"
	end	
	msg = ""
	if nhalf > 0
		msg *= string(nhalf, " half-called genotypes are set to missing. ")
	end
	if nhetero > 0		
		msg *= string(nhetero, " heterozygous genotypes are set to missing")
	end
	if !isempty(msg) 
		verbose && @warn msg
		printconsole(io, verbose, "Warning: "*msg)
	end
	savemagicgeno(magicgenofile,magicgeno)
end

function get_magicped_maxsubnf(magicped::MagicPed; isfounderinbred, model, magicprior)	
	ccpops = MagicBase.get_connected_pops(magicped)
    isempty(ccpops) && error("unexpected connected pops=",ccpops)
	chrmagicprior = magicprior.autosome # despirte of autosome or sex-chromosome
	popmakeup = MagicReconstruct.calpopmakeup(magicped,model, chrmagicprior; isfounderinbred,isautosome=true)
	nfls = [length(i["founder"]) for i=values(popmakeup)]
    maxnf = maximum(nfls)
	maxnf
end

function get_maigcped_connetedcomponents(magicped::MagicPed; isfounderinbred, model, magicprior)	
	ccpops = MagicBase.get_connected_pops(magicped)
    isempty(ccpops) && error("unexpected connected pops=",ccpops)
	chrmagicprior = magicprior.autosome # despirte of autosome or sex-chromosome
	popmakeup = MagicReconstruct.calpopmakeup(magicped,model, chrmagicprior; isfounderinbred,isautosome=true)
	[begin 
		subpops = ccpops[i]
		founders = sort(unique(reduce(vcat, [popmakeup[i]["founder"] for i in subpops])))
		offspring = sort(reduce(vcat, [popmakeup[i]["offspring"] for i in subpops]))
		Dict("subpopulations"=>subpops, "founders"=>founders,"offspring"=>offspring)
	end	for i in eachindex(ccpops)]
end

function aligndiff_fhaplo(bestfhaplo::AbstractDataFrame,runfhaplo::AbstractDataFrame,founders::AbstractVector; isfounderinbred)
	dict = Dict(bestfhaplo[!,:marker] .=> 1:size(bestfhaplo,1))
	posls = [haskey(dict, i) ? dict[i] : nothing for i in runfhaplo[!,:marker]]
	b = .!(isnothing.(posls))
	any(b) && keepat!(posls,b)
	colls = founders .+ 1 # first column is :marker,founders is a list of founder index
	runfhaplo2 = Matrix(runfhaplo[b,colls]) 
	bestfhaplo2 = Matrix(bestfhaplo[posls,colls])
	if isfounderinbred
		perm = MagicBase.permfounder(runfhaplo2,bestfhaplo2)	
		if perm != 1:length(perm)
			runfhaplo2 .= runfhaplo2[:,perm]
		end	
	else		
		# haplotype based permutation
		# runfhaplo3 = phased2haplo(runfhaplo2)
		# bestfhaplo3 = phased2haplo(bestfhaplo2)
		# perm = MagicBase.permfounder(runfhaplo3,bestfhaplo3)	
		# if perm != 1:length(perm)
		# 	runfhaplo3 .= runfhaplo3[:,perm]
		# end	
		# runfhaplo2 = haplo2phased(runfhaplo3)
		# bestfhaplo2 = haplo2phased(bestfhaplo3)
		
		# founder based permutation
		# the absolut phase for each founder (column) of runfhaplo2 is aadjuct to be consistent with bestfhaplo2
		perm = MagicBase.permfounder(runfhaplo2,bestfhaplo2)	
		if perm != 1:length(perm)
			runfhaplo2 .= runfhaplo2[:,perm]
		end	
		MagicBase.alignfounder_phase!(runfhaplo2,bestfhaplo2)		
	end
	posls, bestfhaplo2, runfhaplo2	
end

function phased2haplo(phasedgeno::AbstractMatrix)
    permutedims(reduce(vcat,[reduce(hcat,i) for i in eachcol(phasedgeno)]))
end

function haplo2phased(haplogeno::AbstractMatrix)
	nhaplo = size(haplogeno,2)
	iseven(nhaplo) || @error string("haplo2phased: number of haplotypes is not even; nhaplo=",nhaplo)
    Matrix{Any}(reduce(hcat,[collect(eachrow(haplogeno[:,[i,i+1]])) for i in 1:2:nhaplo-1]))
end

isdiff_allele(x,y) = in("N",[x,y]) ? false : x != y 

isdiff_phasedgeno(x,y) = any(map(isdiff_allele, x,y)) 

isdiff_unphasedgeno(x,y) = any(map(isdiff_allele, x,y)) && any(map(isdiff_allele, reverse(x),y))

function calphasediff_fhaplo(bestfhaplo::AbstractDataFrame,runfhaplo::AbstractDataFrame,founders::AbstractVector; isfounderinbred=false)
	_, bestfhaplo2, runfhaplo2 = aligndiff_fhaplo(bestfhaplo,runfhaplo,founders; isfounderinbred)
	isfounderinbred && error("calphasediff_fhaplo does not work for inbred founders")
	nf = size(bestfhaplo2,2)
	nswitcherr,nswitch = sum(first(MagicBase.calphaseacc(bestfhaplo2[:,i],runfhaplo2[:,i])) for i in 1:nf)
	# @info "calphasediff_fhaplo: phasediff=" (nswitcherr,nswitch)
	nswitch <= 0 ? 0.0 : nswitcherr/nswitch	
end

function caldiff_fhaplo(bestfhaplo::AbstractDataFrame,runfhaplo::AbstractDataFrame,founders::AbstractVector; isfounderinbred)
	_, bestfhaplo2, runfhaplo2 = aligndiff_fhaplo(bestfhaplo,runfhaplo,founders; isfounderinbred)
	if isfounderinbred
		b = map(isdiff_allele, bestfhaplo2, runfhaplo2)
	else
		b = map(isdiff_phasedgeno, bestfhaplo2, runfhaplo2)
	end
	mean(b)
end

function isdiff_fhaplo(bestfhaplo::AbstractDataFrame,runfhaplo::AbstractDataFrame,founders::AbstractVector; isfounderinbred)
	posls, bestfhaplo2, runfhaplo2 = aligndiff_fhaplo(bestfhaplo,runfhaplo,founders; isfounderinbred)
	isdiff = falses(size(bestfhaplo,1), size(bestfhaplo,2)-1)  # first column is :marker
	isdiff2 = view(isdiff,posls, founders)
	if isfounderinbred
		# b = map((x,y)->isdiff_allele(x,y) || (x != "N" && y == "N"), bestfhaplo2, runfhaplo2)
		b = map(isdiff_allele, bestfhaplo2, runfhaplo2)
	else
		b = map(isdiff_phasedgeno, bestfhaplo2, runfhaplo2)
	end
	isdiff2[b] .= true
	isdiff	
end


function impute_refine_chr!(magicgenofile::AbstractString; 
	magicprior::NamedTuple,
    model::AbstractString="jointmodel",
	israndallele::Bool,
	isfounderinbred::Bool=true,	        
	byfounder::Integer, 
	startbyhalf::Integer,
	iscorrectfounder::Bool = true,	
	isimputefounder::Union{Nothing,Bool}=nothing, 
	isallowmissing::Bool,
	threshproposal::Real, 
    isdelmarker::Bool = true,
    delsiglevel::Real = 0.01,
	isinferjunc::Bool = false, 
	isinfererror::Bool = false,	    
	likeparameters::LikeParameters, 
	threshlikeparameters::ThreshLikeParameters, 
	priorlikeparameters::PriorLikeParameters,
	tukeyfence::Real=3.0,		
	minoutlier::Real=0.05, 			
	isspacemarker::Bool = false,
    trimcm::Real=20,
	trimfraction::Real=0.05,  #cM
	skeletonsize::Union{Nothing,Integer},
	isordermarker::Bool = false,
	inputneighbor::Union{Nothing,AbstractDict}=nothing,
	slidewin::Union{Nothing,Integer} = nothing,  
	slidewin_neighbor::Union{Nothing,Integer} = 200,      
	orderactions::AbstractVector = ["inverse","permute"],  
    orderactions_neighbor::AbstractVector = ["inverse11","inverse01"],  
	inittemperature::Real= isordermarker ? 2.0 : 0.0,
    coolrate::Real=0.85,
    minaccept::Real=0.15,
	spacebyviterbi::Bool=false, 
    logio::Union{Nothing,AbstractString,IO},
    imputetempfile::AbstractString,	
	repeatrun::Union{Nothing,Integer}, 	
	maxiter::Integer = 50,	
    verbose::Bool,	
	more_verbose::Bool)
	io = isa(logio,AbstractString) ? open(logio,"w") : logio	
	try	
		filesegs = split(first(MagicBase.split_allext(basename(magicgenofile))),"_")
		tsleepstr = only(filesegs[occursin.("tsleep",filesegs)])
		tsleep = parse(Int,replace(tsleepstr,"tsleep"=>""))
		sleep(tsleep)		
		tused = @elapsed magicgeno = readmagicgeno(magicgenofile)
		chr = 1
	    chrid = magicgeno.markermap[chr][1,:linkagegroup]		
		mem1 = round(Int, memoryuse()/10^6)
		GC.gc()
		mem2 = round(Int, memoryuse()/10^6)				
		msg_run = isnothing(repeatrun) ? "" : string(", run=",repeatrun)
		nmarker = size(magicgeno.markermap[chr],1)
		printconsole(io,verbose,string("chr=", chrid, msg_run, ", reading magicgeno, #markers=", nmarker, ", tused=", round(tused,digits=1),"s",
			", mem=",join([mem1,mem2],"|")))						
		loglikels,fmissls,fmissls_after = impute_refine_chr!(magicgeno,chr;		
			magicprior, model, isfounderinbred, 
			israndallele, inputneighbor, byfounder,	startbyhalf,		
			likeparameters, threshlikeparameters, priorlikeparameters, 
			isinferjunc, iscorrectfounder, isimputefounder, isallowmissing, threshproposal,
			isinfererror, isdelmarker, isspacemarker, isordermarker, 
			tukeyfence, minoutlier, 
			inittemperature, coolrate, delsiglevel, trimcm, trimfraction,
			skeletonsize, slidewin, slidewin_neighbor,
			minaccept,  orderactions,orderactions_neighbor, spacebyviterbi, 
			magicgenofile, imputetempfile,repeatrun, maxiter, verbose, more_verbose,io
		)		
		savemagicgeno(magicgenofile, magicgeno)	 	
		fhaplodf = get_chr_fhaplodf(magicgeno,chr)			
		loglikels, fmissls, fhaplodf
	finally
		isa(logio,AbstractString) && close(io)		
	end
end

function get_chr_fhaplodf(magicgeno::MagicGeno,chr::Integer)
	chrfhaplo = copy(magicgeno.foundergeno[chr])
	founders = magicgeno.magicped.founderinfo[:,:individual]
	fhaplodf = DataFrame(chrfhaplo,founders)
	chrmarkers = copy(magicgeno.markermap[chr][!,:marker])
	insertcols!(fhaplodf, 1, "marker"=>chrmarkers)
	fhaplodf
end

function impute_refine_chr!(magicgeno::MagicGeno,chr::Integer;	    
	magicprior::NamedTuple,
    model::AbstractString, 
	israndallele::Bool,
	isfounderinbred::Bool=true,	    	
    byfounder::Integer,     
	startbyhalf::Integer,
	likeparameters::LikeParameters,
	threshlikeparameters::ThreshLikeParameters,
	priorlikeparameters::PriorLikeParameters,	
	isinferjunc::Bool,	
	iscorrectfounder::Bool, 
	isimputefounder::Union{Nothing,Bool}=nothing, 
	isallowmissing::Bool,
	threshproposal::Real,
	isinfererror::Bool, 
	isdelmarker::Bool, 
	isspacemarker::Bool,
	isordermarker::Bool,
	tukeyfence::Real,		
	minoutlier::Real, 
    delsiglevel::Real,    
    trimcm::Real,
	trimfraction::Real,
    inittemperature::Real,
	coolrate::Real,
	skeletonsize::Union{Nothing,Integer},	
	inputneighbor::Union{Nothing,AbstractDict}=nothing,
	slidewin::Union{Nothing,Integer} = nothing,  
	slidewin_neighbor::Union{Nothing,Integer} = 200,      
	minaccept::Real,		
	orderactions::AbstractVector,
	orderactions_neighbor::AbstractVector,            
	spacebyviterbi::Bool,	
	magicgenofile::AbstractString,
    imputetempfile::AbstractString,
	repeatrun::Union{Nothing,Integer}=nothing, 	
	maxiter::Integer = 50,	
	io::Union{Nothing,IO}, 
    verbose::Bool,	
	more_verbose::Bool=false)	
	popmakeup,priorprocess = MagicReconstruct.calpriorprocess(magicgeno,chr,
		model,magicprior; isfounderinbred)			
	chrid = magicgeno.markermap[chr][1,:linkagegroup]
	issnpGT = [occursin("GT",i) for i in magicgeno.markermap[chr][!,:offspringformat]]
	issnpAD = [occursin("AD",i) for i in magicgeno.markermap[chr][!,:offspringformat]]
	iscalloffspring = model == "depmodel" && any(issnpAD)
	chroffgeno = magicgeno.offspringgeno[chr]
	if iscalloffspring									
		suboffgeno = view(chroffgeno, issnpAD,:)
		callthreshold = 0.95						
		seqerror = MagicBase.get_seqerror(likeparameters)
		for i in eachindex(suboffgeno)
			suboffgeno[i] = MagicBase.callfromprob(MagicBase.genoprobhaplo(suboffgeno[i], seqerror),
				callthreshold; ishaplo=false)
		end
		issnpAD .= false
		issnpGT .= true			
	end	
	mapexpansion = isinferjunc ? magicgeno.magicped.designinfo.mapexpansion : nothing	
	liketargetls, likeerrortuple = get_likeerrortuple(magicgeno,chr; repeatrun, likeparameters,isinfererror,io,verbose)	
	if isnothing(inputneighbor)
		chrneighbor = nothing						
	else
		chrsnpid = magicgeno.markermap[chr][!, :marker]
		neighbor = [get(inputneighbor,i,Vector{String}()) for i in chrsnpid]
		snprule = Dict(chrsnpid .=> 1:length(chrsnpid))
		neighbor2 = [[get(snprule,i,missing) for i in j] for j in neighbor]
		chrneighbor = Dict([i => collect(skipmissing(neighbor2[i])) for i in eachindex(neighbor2)])			
	end		
	# calculate fhaplosetpp
	fhaplosetpp = calfhaploprior(magicgeno,chr)		
	founderformat = unique(magicgeno.markermap[chr][!,:founderformat])	
	isfounderphased = issubset(["GT_phased"], founderformat)

	# modfied: priorprocess
	loglikels, chrfhaplo, snporder, offexcl, likeerrortuple,correctdf,fmissls,fmissls_after = impute_refine_chr!(magicgeno.magicped, chroffgeno, chrneighbor,
		popmakeup, priorprocess, fhaplosetpp;
		isfounderinbred, isfounderphased, byfounder, startbyhalf,
		israndallele, issnpGT,issnpAD,  chrid,  			
		iscorrectfounder, isimputefounder, isallowmissing, threshproposal, isinfererror, isdelmarker, isspacemarker,isordermarker,
		liketargetls, likeerrortuple, threshlikeparameters, priorlikeparameters, tukeyfence,minoutlier, 
		inittemperature, coolrate, delsiglevel, trimcm, trimfraction,
		skeletonsize, slidewin, slidewin_neighbor,
		minaccept,  orderactions,orderactions_neighbor, 
		isinferjunc, mapexpansion, spacebyviterbi, imputetempfile,repeatrun, 
		maxiter, verbose, more_verbose,io
	)		
	# @info correctdf
	if  !isempty(correctdf)
		correctdf2 = transferchangedf(correctdf,magicgeno.markermap[chr],
			magicgeno.magicped.founderinfo)
		MagicBase.pushmisc!(magicgeno,"correction"=>correctdf2)
	end
	save_impute_refine_chr!(magicgeno, chr, chrfhaplo,likeerrortuple,offexcl; isfounderinbred,isinfererror)
	MagicReconstruct.check_prior_consistency(magicgeno,chr, priorprocess,snporder)
	if iscalloffspring					
		oldmagicgeno = readmagicgeno(magicgenofile)
		magicgeno.offspringgeno[chr] = oldmagicgeno.offspringgeno[chr]			
	end
	MagicReconstruct.setmagicgeno!(magicgeno, chr, priorprocess,snporder)  # delete and order markers 	
	loglikels,fmissls,fmissls_after
end



function impute_refine_chr!(magicped::MagicPed, chroffgeno::AbstractMatrix,
	chrneighbor::Union{Nothing,AbstractDict},
    popmakeup::AbstractDict, priorprocess::AbstractDict,
    fhaplosetpp::AbstractVector;    	
	isfounderinbred::Bool,
	isfounderphased::Bool, 
	byfounder::Integer, 
	startbyhalf::Integer,
	chrid::AbstractString,         	    
    israndallele::Bool=true,     
	issnpGT::AbstractVector, 
	issnpAD::AbstractVector, 
	liketargetls::AbstractVector,
	likeerrortuple::NamedTuple, 
	threshlikeparameters::ThreshLikeParameters,
	priorlikeparameters::PriorLikeParameters,	
	isimputefounder::Union{Nothing, Bool},
	isallowmissing::Bool, 
	threshproposal::Real,
	iscorrectfounder::Bool, 
	isinfererror::Bool, 
	isdelmarker::Bool, 
	isspacemarker::Bool=false,
	isordermarker::Bool=false, 
	tukeyfence::Real=3.0,		
	minoutlier::Real=0.05, 
    delsiglevel::Real=0.01,    
	trimcm::Real=20,
	trimfraction::Real=0.05, 
    inittemperature::Real=0.0,
	coolrate::Real=0.85,
	skeletonsize::Union{Nothing,Integer},		
	slidewin::Union{Nothing,Integer} = nothing,  
	slidewin_neighbor::Union{Nothing,Integer} = 200,      
	minaccept::Real=0.1,		
	orderactions::AbstractVector=[],
	orderactions_neighbor::AbstractVector=[],                
	isinferjunc::Bool=false,
	mapexpansion::Union{Nothing,Real}=nothing, 
	spacebyviterbi::Bool, 
    imputetempfile::AbstractString,
	repeatrun::Union{Nothing,Integer}=nothing, 		
	maxiter::Integer = 50,	
	io::Union{Nothing,IO}=nothing, 
    verbose::Bool=true,	
	more_verbose::Bool=false)	
	chrstartt = time()
	# founder partition and initialization	
	minfmiss = isallowmissing ? 1e-4 : 0.0
	chrfhaplo,fmissls, isimputefounder = MagicImpute.get_initfhaplo(fhaplosetpp; minfmiss, isimputefounder)	    
	avgfmiss = mean(fmissls)
	defaultby = MagicBase.get_partition_defaultby(isfounderinbred,isfounderphased)	
	findexlist = MagicBase.getfindexlist(byfounder,fmissls, popmakeup;defaultby,minfmiss)			
	msg_repeatrun = isnothing(repeatrun) ? "" : string(", run=", repeatrun)
	msg = string("chr=", chrid, msg_repeatrun)
	if isimputefounder		
		# avgfmissblock = [mean(fmissls[findex]) for findex in findexlist]					
		msg *= string(", partition=", join.(findexlist,"|"))		
		# 	", fmiss_blocks=", round.(avgfmissblock,digits=4))
		msg *= string(", avgfmiss=", round(avgfmiss,digits=4), ", fmiss=", round.(fmissls,digits=4))		
	else
		msg *= string(", avgfmiss=", round(avgfmiss,digits=4),	", fmiss=", round.(fmissls,digits=4), ", skip imputing founders")
	end	
	printconsole(io,verbose,msg)		
	
	# initialization
	nsnp = size(chroffgeno,1)
	ismalexls = MagicImpute.getismalexls(magicped,chrid)
	priorspace= MagicReconstruct.getpriorstatespace(magicped; isfounderinbred)				
	founder2progeny = index_founder2progeny(magicped)	
	# maxwinsize = div(nsnp,min(20,5+round(Int, nsnp/500)))
	if isnothing(slidewin_neighbor) 
		if nsnp < 500			
			maxwinsize = max(10,min(50,div(nsnp,5)))
		elseif nsnp < 2000
			maxwinsize = div(nsnp, 10)
		else
			maxwinsize = 200
		end
	else
		maxwinsize = slidewin_neighbor  # default 200
		maxwinsize = max(10,min(div(nsnp,5),maxwinsize))
	end		
	nbrmaxwin = [maxwinsize for _ in 1:4]
	if isnothing(slidewin)
		if isnothing(chrneighbor)				
			initwinsize = max(10,round(Int, sqrt(nsnp)))
		else
			chrnneighbor = mean(length.(values(chrneighbor)))
			initwinsize = max(10,round(Int,min(sqrt(nsnp),chrnneighbor/3)))				
		end
	else
		initwinsize = slidewin
	end							
	initwinsize = min(maxwinsize,max(2, initwinsize))
	#  set winsizels for ordering marker
	if inittemperature <= 0.01
		niter = 4
	else
		niter = inittemperature <= 0.5 ? 0 : ceil(Int,-1*log(inittemperature/0.5)/log(coolrate))+1 # interstions from T0 to 0.5
		niter += 6 # 6 iterations from (<0.5) to 0
	end
	winsizels = get_winsizels(initwinsize,niter)
	last(winsizels) <= 3 || @error string("unexpected winsizels=",winsizels)
	pushfirst!(winsizels,first(winsizels)) # first iteration does not perform ordering
	miniteration_order = min(maxiter, length(winsizels)+1)
	if length(winsizels) < maxiter
		append!(winsizels,repeat([last(winsizels)], maxiter-length(winsizels)))
	end	
	offspringexcl = []	    
	# offlogl = -Inf*ones(size(chroffgeno,2))
	# priorlength = min(20.0, max(5.0,round(2000/nsnp,digits=1)))
	priorlength = min(10.0, round(200/nsnp,digits=1))
	chrlenls = [-Inf, -Inf]
	ndiffls = Vector{Int}()
	ncorrectls = Vector{Int}()	  
	merrorlsls = Vector{Vector{Float64}}()
	msgwin = ""
	if isordermarker 
		if !isnothing(chrneighbor)				
			msgwin *= string(", max_slidewin_neighbor=",maxwinsize) 				
		end	
		msgwin *= string(", init_slidewin=",initwinsize)					
	end
	msg_repeatrun = isnothing(repeatrun) ? "" : string(", run=", repeatrun)
	msg = string("chr=", chrid, msg_repeatrun, ", #marker=", nsnp, msgwin, ", begin")	
	printconsole(io,verbose,msg)		
	snporder = collect(1:nsnp)		
	temperature = inittemperature	
	step_verbose = false					
	isinferseq = isinfererror && any(issnpAD)			    
	actionlsls = [[isimputefounder, iscorrectfounder, isinfererror, isinferseq, isdelmarker, isspacemarker,isordermarker]]
	oldmapexpansion = mapexpansion			
	correctdf = DataFrame()	
	upbyhalf = false
	alwaysacceptls = [true,true]
	imputestuckls = [-1,-1]
	for it in 1:maxiter						
		slidewinsize = truncated(Poisson(max(2, winsizels[it])),2,maxwinsize)
		reversechr = iseven(it)
		# re-randomize findexlist reults in worse imputation for large byfounder and large #missingparnts
		# if first(last(actionlsls))	
		# 	# byfounder2 = last(ndiffls)/length(chrfhaplo) < 0.05 ? ceil(Int, byfounder/2) : byfounder			
		# 	findexlist = MagicBase.getfindexlist(byfounder,fmissls, popmakeup;defaultby,minfmiss)			
		# end
		# modify chrfhaplo, priorprocess,chrlenls, ndiffls,ncorrectls, merrorlsls, actionlsls, offspringexcl
		likeerrortuple, offspringexcl, correctdf0, tused, msg, logbook_order, upbyhalf = impute_refine_chr_it!(chrfhaplo,chroffgeno,
			popmakeup,priorprocess, priorspace, fhaplosetpp;
			isfounderinbred, israndallele, ismalexls, founder2progeny,findexlist, 
			likeerrortuple, snporder, issnpGT, 
			threshlikeparameters, priorlikeparameters, liketargetls, tukeyfence,minoutlier, offspringexcl, 											
			temperature, reversechr,isallowmissing, threshproposal, startbyhalf, 
			delsiglevel, priorlength, trimcm, trimfraction, 
			chrneighbor, minaccept,nbrmaxwin, slidewinsize, maxwinsize, orderactions,orderactions_neighbor, 
			chrlenls, ndiffls, ncorrectls, merrorlsls, actionlsls, 
			alwaysacceptls, imputestuckls, upbyhalf, 
			iteration=it, miniteration_order, imputetempfile, spacebyviterbi, step_verbose)			
		if isinferjunc
			startt_junc = time()
			juncdist = infer_juncdist(chrfhaplo,chroffgeno,popmakeup,priorprocess,priorspace;
				likeerrortuple, decodetempfile = imputetempfile, israndallele, issnpGT,snporder);
			magicgeno.magicped.designinfo = juncdist
			MagicBase.reset_juncdist!(magicgeno.magicped,model; verbose=false,isfounderinbred,isinferjunc)
			magicprior=MagicReconstruct.calmagicprior(magicgeno,model;isfounderinbred);
			popmakeup = MagicReconstruct.calpopmakeup(magicgeno,chr, model,magicprior; isfounderinbred)
			priorprocess = MagicReconstruct.uppriorprocess(priorprocess,popmakeup)
			(;mapexpansion) = juncdist
			msg *= string(", R=",mapexpansion)				
			push!(tused,string(round(Int,time()-startt_junc)))
			isinferjunc = abs(mapexpansion  - oldmapexpansion) > 0.001
			oldmapexpansion = mapexpansion				
		end
		if !isempty(correctdf0)
			correctdf = vcat(correctdf,correctdf0)
		end		
		msg = string("chr=", chrid, msg_repeatrun, ", it=",it,msg)			
		inittemperature > 0 && (msg *= string(", T=", round(temperature,digits=3)))			
		reversechr = iseven(it+1)	
		msg *= string(", t=", join(tused,"|"),"s")
		# mem1 = round(Int, memoryuse()/10^6)
		# GC.gc()
		# mem2 = round(Int, memoryuse()/10^6)
		# msg *= string(", mem=",mem1,"|",mem2, "MB")	        
		if isempty(tused) 				
			if any(last(actionlsls)) && it < maxiter  
				msg *= string(",unexpected actionls=",last(actionlsls))
				@warn msg
				printconsole(io,false,"Warning: "*msg)
			end
		else
			printconsole(io,verbose,msg)
		end
		if !isnothing(logbook_order)
			if step_verbose
				show(stdout, "text/plain", logbook_order)
				print("\n")
			end
			if more_verbose
				show(io, "text/plain", logbook_order)
				write(io, "\n")
			end
			# first row of logbook_order: ["proptype", "win","nprop","tused","apt","logl"]
			israndorder = occursin.(r"_Poisson$",logbook_order[:,1])
			if any(israndorder)
				avgaccept = mean(logbook_order[israndorder,5])
				if avgaccept < minaccept						
					push!(winsizels, last(winsizels))
					popfirst!(winsizels)
				end
			end
		end		
		if it >=2				
			if temperature<=1e-4
				temperature = 0.0						
			elseif temperature<=0.1
				temperature *= 0.1			
			else
				temperature *= temperature <= 0.5 ? max(0.5,coolrate^3) : coolrate					
			end
		end
		if it == maxiter
			msg = string("chrid=", chrid, ",it=", it, ", reach the max iteration!")
			@warn msg
			printconsole(io, false, msg)
		end										
		if !any(last(actionlsls)) || it == maxiter				
			if isspacemarker
				# prior2df(prior) = DataFrame(deltd = first(values(prior)).markerdeltd,incl = first(values(prior)).markerincl)					
				# befpriorprocess = deepcopy(priorprocess)
				pri1 = first(values(priorprocess))					
				nsnpincl = sum(pri1.markerincl)	
				ndiffloc = sum([!isnothing(i) && i > 0 for i in pri1.markerdeltd])				
				isskeleton = isnothing(skeletonsize) ? ndiffloc > 20 : (nsnpincl > skeletonsize)
				if isskeleton						
					skeletonprior = markerskeleton_chr(chrfhaplo,chroffgeno, popmakeup,priorprocess;
						israndallele, isinfererror, priorlikeparameters,liketargetls, 
						likeerrortuple, offspringexcl, 							
						snporder,issnpGT, isinferseq, 
						decodetempfile = imputetempfile, skeletonsize,  
						chrid, startit=it+1, spacebyviterbi, verbose, io)
					rescalemap!(priorprocess,skeletonprior)						
				else
					msg  = string("chr=", chrid, msg_repeatrun, ", #marker=", nsnp, ", #marker_incl=", nsnpincl,
						", #segment=",ndiffloc, ", no distance rescale by skeleton markers")
					printconsole(io,verbose, msg)
				end
			end
			break
		end			
	end				
	loglikels = cal_hmm_loglikels(chrfhaplo,chroffgeno,popmakeup,priorprocess;
		likeerrortuple, offspringexcl, decodetempfile=imputetempfile, 
		israndallele,issnpGT,snporder, incldelmarkers = true, 
	)
	ncorrect = size(correctdf,1)
	pri1=first(values(priorprocess))
	snpincl = snporder[pri1.markerincl]			
	nmissallele = sum(chrfhaplo[snpincl, :] .== "N")	
	# !isfounderinbred || (nmissallele = div(nmissallele, 2))
	fmissls_after = mean(chrfhaplo[snpincl,:] .== "N",dims=1)[1,:]
	if !isfounderinbred
		fmissls_after = mean.(Iterators.partition(fmissls_after,2))	
	end
	nsnpincl = sum(pri1.markerincl)		
	chrlen = 100*sum(pri1.markerdeltd[pri1.markerincl][1:end-1])
	mem1 = round(Int, memoryuse()/10^6)
	GC.gc()
	mem2 = round(Int, memoryuse()/10^6)		
	msg  = string("chr=", chrid, msg_repeatrun, 
		", logl=", round(sum(loglikels),digits=1),		
		", #marker=", nsnp, 
		", #marker_incl=", nsnpincl,
		", #missf=", nmissallele,
		iscorrectfounder ? string(", #correctf=", ncorrect) : "",			
		", l=", round(chrlen,digits=1), "cM",		
		", fmiss=", join(round.(fmissls_after,digits=4),"|"), 	
		", t=", round(Int, time()-chrstartt), "s",			
		", mem=",mem1,"|",mem2, "MB", 
		", end", 
	)
	printconsole(io,verbose,msg)		
	loglikels, chrfhaplo, snporder, offspringexcl, likeerrortuple, correctdf, fmissls,fmissls_after
end	

function cal_hmm_loglikels(chrfhaplo::AbstractMatrix,chroffgeno::AbstractMatrix,        
    popmakeup::AbstractDict,priorprocess::AbstractDict;    
    likeerrortuple::NamedTuple,     
    snporder::AbstractVector,    
    decodetempfile::AbstractString,    
    israndallele::Bool,
	incldelmarkers::Bool, 
    offspringexcl::AbstractVector, 
    issnpGT::AbstractVector)
    (;epsfls, epsols, epsols_perind, seqerrorls,allelebalancemeanls,allelebalancedispersels, alleledropoutls) = likeerrortuple	
    loglikels = MagicReconstruct.hmm_loglikels(chrfhaplo,chroffgeno,popmakeup,priorprocess;
        epsf=epsfls, epso=epsols,  epso_perind = epsols_perind, 
        seqerror = seqerrorls, allelebalancemean = allelebalancemeanls, 
        allelebalancedisperse = allelebalancedispersels, alleledropout = alleledropoutls, 
        decodetempfile, israndallele,issnpGT,snporder    
	)	
	if incldelmarkers
		delsnps = snporder[.!(first(values(priorprocess)).markerincl)]		    
		if !isempty(delsnps) 
			singlelogl = MagicImpute.calsinglelogl(chrfhaplo,chroffgeno;
				snpincl = delsnps, popmakeup,
				epsf=epsfls, epso=epsols,  epso_perind = epsols_perind, 
				seqerror = seqerrorls, allelebalancemean = allelebalancemeanls, 
				allelebalancedisperse = allelebalancedispersels, alleledropout = alleledropoutls, 
				israndallele,issnpGT
			)
			loglikels .+= sum(singlelogl,dims=1)[1,:]
		end
	end
	loglikels[offspringexcl] .= 0.0 # loglike is set to 0 for excluded offspring
    loglikels
end

# function set_fhaplosetpp_missing!(fhaplosetpp::AbstractVector, chrfhaplo::AbstractMatrix, bdiff::AbstractMatrix, snpincl::AbstractVector)
# 	isfounderinbred = length(fhaplosetpp) == size(bdiff,2)	
# 	size(chrfhaplo, 1) == size(bdiff,1) == length(snpincl) || @error "inconsistent #markers"
# 	if isfounderinbred
# 		for i in findall(bdiff)
# 			(snp0, p) = Tuple(i)		
# 			snp = snpincl[snp0]	
# 			fhaplosetpp[p][snp] = ["1","2"]
# 		end
# 	else
# 		bdidff_outbred = permutedims(reduce(hcat,[Vector.(collect(Iterators.partition(i,2))) for i in eachrow(bdiff)]))
# 		for i in eachindex(IndexCartesian(),bdidff_outbred)
# 			bdidff_outbred[i] == [false,false] && continue
# 			(snp0, p) = Tuple(i)			
# 			snp = snpincl[snp0]
# 			if bdidff_outbred[i] == [true,true] 							
# 				fhaplosetpp[p][snp] = [["1","1"],["1","2"],["2","1"],["2","2"]]
# 			elseif bdidff_outbred[i] == [true,false] 							
# 				a = chrfhaplo[snp0,2*p] # use snp0 inst4ead of snp, because chrfhaplo and bdiff correspoind to the same markers (snpincl)
# 				fhaplosetpp[p][snp] = union(fhaplosetpp[p][snp], [["1",a],["2",a]])
# 			elseif bdidff_outbred[i] == [false,true] 
# 				a = chrfhaplo[snp0,2*p-1]
# 				fhaplosetpp[p][snp] = union(fhaplosetpp[p][snp], [[a,"1"],[a,"2"]])
# 			end
# 		end
# 	end
#     fhaplosetpp
# end

# modify chrfhaplo, priorprocess,chrlenls, ndiffls, actionlsls
# modified fhaplosetpp if iscorrectfounder =true
function impute_refine_chr_it!(chrfhaplo::AbstractMatrix, chroffgeno::AbstractMatrix,
    popmakeup::AbstractDict,
    priorprocess::AbstractDict, priorspace::AbstractDict,
    fhaplosetpp::AbstractVector;    	
	isfounderinbred::Bool, 
    ismalexls::AbstractVector,
    founder2progeny::AbstractVector,
    findexlist::AbstractVector, 	
	liketargetls::AbstractVector,
	likeerrortuple::NamedTuple,     
    threshlikeparameters::ThreshLikeParameters,
	priorlikeparameters::PriorLikeParameters,	
	offspringexcl::AbstractVector,
	snporder::AbstractVector,
    israndallele::Bool=true,     
	issnpGT::AbstractVector, 	
	tukeyfence::Real,	
	minoutlier::Real,			
    reversechr::Bool,
    delsiglevel::Real,
    priorlength::Real,
    trimcm::Real,
	trimfraction::Real, 
    temperature::Real,
	chrneighbor::Union{Nothing, AbstractDict},
    slidewinsize::Union{Integer,Distribution},
	maxwinsize::Integer,
	minaccept::Real,	
	nbrmaxwin::AbstractVector,
	orderactions::AbstractVector,
	orderactions_neighbor::AbstractVector,
    chrlenls::AbstractVector,
    ndiffls::AbstractVector,
    ncorrectls::AbstractVector,
    merrorlsls::AbstractVector, 
	actionlsls::AbstractVector,		
	imputestuckls::AbstractVector, 
	alwaysacceptls::AbstractVector, 
	upbyhalf::Bool, 	
	isallowmissing::Bool,
	threshproposal::Real,
	iteration::Integer,	
	miniteration_order::Integer,
	startbyhalf::Integer, 
    imputetempfile::AbstractString,
	spacebyviterbi::Bool,
	step_verbose::Bool=false)
	# deletion after infereps step since small inital eps results in too many deletions
	# correction after 2 iterations of deletion since  correction may change mis-groupped markers.	
    isimputefounder, iscorrectfounder, isinfererror, isinferseq, isdelmarker, isspacemarker, isordermarker = last(actionlsls)	
	(;epsfls, epsols, epsols_perind, seqerrorls,allelebalancemeanls,allelebalancedispersels, alleledropoutls) = likeerrortuple
    tused =[]
    msg = ""		
	errortuples = (epsf=epsfls,epso=epsols,  epso_perind = epsols_perind, seqerror = seqerrorls, allelebalancemean=allelebalancemeanls, 
		allelebalancedisperse=allelebalancedispersels, alleledropout = alleledropoutls)		
	miditeration = max(5, startbyhalf)
    if isimputefounder 
        startt = time()			
		alwaysaccept = alwaysacceptls[end]
		imputestuck = imputestuckls[end]
		if upbyhalf	&& imputestuck >= 3
			isimputefounder = false							
		end	
		if !isallowmissing && !isimputefounder # last iteration to impute all missing founder genotypes			
			threshproposal2 = -1
			alwaysaccept2 = true			
		else
			threshproposal2 = threshproposal
			alwaysaccept2 = alwaysaccept			
		end
        deltloglike, ndiff = founderimpute_chr!(chrfhaplo,chroffgeno, popmakeup,priorprocess, fhaplosetpp;
            findexlist, errortuples..., offspringexcl, snporder, 
			reversechr = iseven(iteration), 
			alwaysaccept = alwaysaccept2, upbyhalf, threshproposal = threshproposal2, 
			israndallele,issnpGT,imputetempfile)							
        push!(ndiffls,ndiff)		
		if alwaysaccept2		
			imputestuck = -1
		else
			if deltloglike > 0 && ndiff > 0
				imputestuck = 0
			else
				deltstuck = (deltloglike <= 0.0) + (ndiff== 0)
				imputestuck =  max(0, imputestuck) + deltstuck				
			end			
		end
		push!(imputestuckls,imputestuck)
        msg *= string(", #diff=",ndiff, ", Δlogl=",round(deltloglike,digits=1))
		msg *= string(", stuck=",imputestuck==-1 ? -Inf : imputestuck, ", byhalf=",upbyhalf)		
		snpincl = snporder[first(values(priorprocess)).markerincl]		
		nmissallele = sum(chrfhaplo[snpincl, :] .== "N")
		# isfounderinbred || (nmissallele = div(nmissallele, 2))
		msg *= string(", #missf=",nmissallele)		
        step_verbose && println(msg)
        push!(tused,string(round(Int,time()-startt)))		

		# update alwaysaccept and upbyhalf
		if alwaysaccept							
			if ndiffls[end] == 0
				alwaysaccept = false
			elseif upbyhalf && iteration >= miditeration
				if  ndiffls[end]/length(chrfhaplo) < 0.001
					alwaysaccept = false							
				elseif length(ndiffls) >=3 && allequal(ndiffls[end-2:end])					
					alwaysaccept = false		
				elseif length(ndiffls) >= 4 && sum(ndiffls .<= ndiffls[end]) >= 4										
					alwaysaccept = false			
				end
			end
			if alwaysaccept
				if !upbyhalf 					
					upbyhalf = iteration >= (startbyhalf - 1)
				end
			else
				upbyhalf = true				
			end						
		else
			upbyhalf = true
		end		
		push!(alwaysacceptls,alwaysaccept)
    end	

	if iscorrectfounder && (iteration >= 2 || !isimputefounder)		
		startt = time()				
		correctdf = foundercorrect_chr!(chrfhaplo, chroffgeno,            
			popmakeup,priorprocess,priorspace, fhaplosetpp;
			errortuples..., offspringexcl, snporder, 
			decodetempfile = imputetempfile,
			inclmissingallele = false,  # if true, most of them will be missing in the following iterations
			ismalexls,founder2progeny,israndallele, issnpGT, step_verbose)
		ncorrect =  isempty(correctdf) ? 0 : size(correctdf,1)
		push!(ncorrectls,ncorrect)		
		if (isimputefounder && upbyhalf) || !isimputefounder
			if ncorrect == 0 
				iscorrectfounder = false        			
			elseif length(ncorrectls) >= 3 && allequal(ncorrectls[end-2:end]) 
				iscorrectfounder = false        			
			elseif  isimputefounder && length(ncorrectls) >= 4 && all(ndiffls[end-2:end] .== ncorrectls[end-3:end-1])
				# imputation and correction are locked with each other
				iscorrectfounder = false        							
			end		
		end
		msg *= string(", #correctf=",ncorrect)
		step_verbose && println(msg)
		push!(tused,string(round(Int,time()-startt)))
	else
		correctdf = DataFrame()
	end	
	if 1 < iteration <= miditeration && isimputefounder
		monosnps = findall([unique(i) in [["1"],["2"]] for i in eachrow(chrfhaplo)])
		if !isempty(monosnps)
			msg *= string(", #mono=",length(monosnps))
			genoset = length(fhaplosetpp) == size(chrfhaplo,2) ? ["1","2"] : [["1", "1"], ["1", "2"], ["2", "1"],["2", "2"]]
			for i in eachindex(fhaplosetpp)
				fhaplosetpp[i][monosnps] .= [genoset for _ in eachindex(monosnps)]			
			end
		end
	end
	if iteration > miditeration || !isimputefounder
        snpincl = snporder[first(values(priorprocess)).markerincl]		
		monosnps = findall([unique(i) in [["1"],["2"]] for i in eachrow(chrfhaplo)])
		intersect!(monosnps,snpincl)
		if !isempty(monosnps) && length(snpincl) - length(monosnps) >= 2 
			msg *= string(", #del_mono=",length(monosnps))
			MagicReconstruct.setpriorprocess!(priorprocess, snporder,monosnps)			
		end		        
	end
	if isinfererror && (iteration >= 2 || !isimputefounder)		
		merrorls = zeros(7)
		absthreshls = [1e-3,1e-3,1e-3,1e-3,1e-2,5e-3,1e-3]
        startt = time()				 
		bool = length(merrorlsls) >= 2 
		targetls = intersect(["foundererror","offspringerror"],liketargetls)		
		pri1=first(values(priorprocess))
        snpincl = snporder[pri1.markerincl]		   		
		infer_err_fun = infer_errorls_viterbi! 
		infer_err_fun(chrfhaplo,chroffgeno, popmakeup,priorprocess;
			targetls, priorlikeparameters,
			epsf=epsfls,epso=epsols,  epso_perind = epsols_perind, 
			seqerror = seqerrorls, allelebalancemean=allelebalancemeanls, allelebalancedisperse=allelebalancedispersels,
			alleledropout = alleledropoutls, offspringexcl, snporder, 
			decodetempfile = imputetempfile, israndallele, issnpGT, temperature, itmax = iteration <= 4 ? 20 : 10)		
		if in("foundererror",targetls)
			errindex = 1
			merrorls[errindex]  = isa(epsfls,AbstractVector) ? mean(epsfls[snpincl]) : epsfls      
			if length(merrorlsls) >=2 							
				bool *= abs(merrorlsls[end-1][errindex] - merrorls[errindex]) <= absthreshls[errindex] || abs(merrorlsls[end][errindex] - merrorls[errindex]) <= absthreshls[errindex]
			end	
			msg *= string(", εf=", round(mepsf,digits=3))			
		end
		errindex = 2
		merrorls[errindex] = isa(epsols,AbstractVector) ? mean(epsols[snpincl]) : epsols     
		if length(merrorlsls) >=2 							
			bool *= abs(merrorlsls[end-1][errindex] - merrorls[errindex]) <= absthreshls[errindex] || abs(merrorlsls[end][errindex] - merrorls[errindex]) <= absthreshls[errindex]
		end	
		msg *= string(", εo=", round(merrorls[errindex],digits=3))					
		if in("peroffspringerror",liketargetls) 				
			infer_error_perind_viterbi!(chrfhaplo,chroffgeno, popmakeup,priorprocess;
				epsf=epsfls,epso=epsols,  epso_perind = epsols_perind, 
				seqerror = seqerrorls, allelebalancemean=allelebalancemeanls, 
				allelebalancedisperse=allelebalancedispersels,alleledropout = alleledropoutls, 
				prior_peroffspringerror = priorlikeparameters.peroffspringerror,         				
				israndallele, issnpGT,snporder,decodetempfile = imputetempfile, 
				temperature, itmax = iteration <= 4 ? 20 : 10)				
			errindex = 3
			merrorls[errindex] = mean(epsols_perind)
			if length(merrorlsls) >=2 							
				bool *= abs(merrorlsls[end-1][errindex] - merrorls[errindex]) <= absthreshls[errindex] || abs(merrorlsls[end][errindex] - merrorls[errindex]) <= absthreshls[errindex]
			end	
			msg *= string(", ξo=", round(merrorls[errindex],digits=3))						
		end
		if isinferseq
			snpADls = intersect(snpincl,findall(.!issnpGT))			
			targetls = setdiff(liketargetls,["foundererror","offspringerror","peroffspringerror"])			
			infer_err_fun(chrfhaplo,chroffgeno, popmakeup,priorprocess;
				targetls, priorlikeparameters,					
				epsf=epsfls, epso=epsols,  epso_perind = epsols_perind, 
				seqerror = seqerrorls,allelebalancemean = allelebalancemeanls, 
				allelebalancedisperse = allelebalancedispersels, alleledropout = alleledropoutls, 
				offspringexcl, israndallele, issnpGT, snporder,decodetempfile = imputetempfile,
				temperature,itmax = iteration <= 4 ? 20 : 10)		
			merror = isa(seqerrorls,AbstractVector) ? mean(seqerrorls[snpADls]) : seqerrorls
			mbalancemean = isa(allelebalancemeanls,AbstractVector) ? mean(allelebalancemeanls[snpADls]) : allelebalancemeanls
			mbalancedisperse = isa(allelebalancedispersels,AbstractVector) ? mean(allelebalancedispersels[snpADls]) : allelebalancedispersels
			mdropout = isa(alleledropoutls,AbstractVector) ? mean(alleledropoutls[snpADls]) : alleledropoutls
			merrorls[4:7] .= [merror,mbalancemean,mbalancedisperse,mdropout]				
			if length(merrorlsls) >=2 											
				for errindex in 4:7
					bool *= abs(merrorlsls[end-1][errindex] - merrorls[errindex]) <= absthreshls[errindex] || abs(merrorlsls[end][errindex] - merrorls[errindex]) <= absthreshls[errindex]
				end
			end			
			in("seqerror",liketargetls) && (msg *= string(", εseq=",round(merrorls[4],digits=4)))
			in("allelebalancemean",liketargetls) &&  (msg *= string(", ABmean=",round(merrorls[5],digits=2)))
			in("allelebalancedisperse",liketargetls) && (msg *= string(", disperse=",round(merrorls[6],digits=3)))
			in("alleledropout",liketargetls) && (msg *= string(", dropout=",round(merrorls[7],digits=3)))
		end		
        if bool && ((iteration > miditeration && temperature <= 0.5) || !isimputefounder)			
			if bool && !any([isimputefounder,iscorrectfounder,isdelmarker, temperature > 0])	
				# bool && !any([isimputefounder,iscorrectfounder,isdelmarker, isordermarker])	
				if all(@. merrorls[1:3] - merrorlsls[end][1:3] <= 1e-3)
					isinfererror = false	
				end
			end
			delsnps_epsf = []	
			delsnps_epso = []	
			delsnps_seqerr = []				
			delsnps_dropout = []
			delsnps_abm = []
			delsnps_abd = []
			if isa(epsfls,AbstractVector) 				
				delsnps_epsf= MagicReconstruct.del_error_indices(epsfls[snpincl];
					tukeyfence,minoutlier=minoutlier,max_error=threshlikeparameters.foundererror)
				delsnps_epsf = snpincl[delsnps_epsf]
			end	
			if isa(epsols,AbstractVector) 				
				delsnps_epso= MagicReconstruct.del_error_indices(epsols[snpincl];
					tukeyfence,minoutlier=minoutlier,max_error=threshlikeparameters.offspringerror)
				delsnps_epso = snpincl[delsnps_epso]
			end		
			if isinferseq && isa(seqerrorls,AbstractVector) && in("seqerror",liketargetls)
				delsnps_seqerr= MagicReconstruct.del_error_indices(seqerrorls[snpADls];
					tukeyfence,minoutlier=minoutlier,max_error=threshlikeparameters.seqerror)
				delsnps_seqerr = snpADls[delsnps_seqerr]
			end				
			if isinferseq && isa(allelebalancemeanls,AbstractVector) && in("allelebalancemean",liketargetls) 
				abmls = allelebalancemeanls[snpADls]					
				max_balancemean = threshlikeparameters.allelebalancemean
				delsnps_abm = findall(@. abmls > max_balancemean || abmls < 1-max_balancemean)
				snps_outlier = MagicReconstruct.get_outerlier_indices(abmls; tukeyfence, side="both")
				union!(delsnps_abm,snps_outlier)
				sort!(delsnps_abm)
				delsnps_abm = snpADls[delsnps_abm]
			end	
			if isinferseq && isa(allelebalancedispersels,AbstractVector) && in("allelebalancedisperse",liketargetls) 
				abdls = allelebalancedispersels[snpADls]					
				delsnps_abd = findall(abdls .> threshlikeparameters.allelebalancedisperse)				
				delsnps_abd = snpADls[delsnps_abd]
			end	
			if isinferseq && isa(alleledropoutls,AbstractVector) && in("alleledropout",liketargetls) 	
				adropls = alleledropoutls[snpADls]
				delsnps_dropout = findall(adropls .> threshlikeparameters.alleledropout)				
				delsnps_dropout = snpADls[delsnps_dropout]
			end	
			delsnps = union(delsnps_epsf,delsnps_epso,delsnps_seqerr,delsnps_abm,delsnps_abd,delsnps_dropout)
			if !isempty(delsnps) && length(snpincl) - length(delsnps) >= 2		
				isinfererror = true		
				MagicReconstruct.setpriorprocess!(priorprocess,snporder, delsnps)
				msg *= string(", #del_err=",length(delsnps))	
				if !isempty(delsnps_epsf)
					msg *= string(", del_εf=",join(round.(epsfls[delsnps_epsf],digits=2),"|"), " at snps=",join(delsnps_epsf, "|"))				
				end
				if !isempty(delsnps_epso)
					msg *= string(", del_εo=",join(round.(epsols[delsnps_epso],digits=2),"|"), " at snps=",join(delsnps_epso, "|"))				
				end
				if !isempty(delsnps_seqerr)
					msg *= string(", del_εseq=",join(round.(seqerrorls[delsnps_seqerr],digits=3),"|"), " at snps=",join(delsnps_seqerr, "|"))				
				end
				if !isempty(delsnps_abm)
					msg *= string(", del_ABm=",join(round.(allelebalancemeanls[delsnps_abm],digits=3),"|"), " at snps=",join(delsnps_abm, "|"))				
				end
				if !isempty(delsnps_abd)
					msg *= string(", del_disperse=",join(round.(allelebalancedispersels[delsnps_abd],digits=3),"|"), " at snps=",join(delsnps_abd, "|"))				
				end	
				if !isempty(delsnps_dropout)
					msg *= string(", del_dropout=",join(round.(alleledropoutls[delsnps_dropout],digits=3),"|"), " at snps=",join(delsnps_dropout, "|"))				
				end
			end
        end
		isinferseq = isinferseq && isinfererror
        step_verbose && println(msg)
        push!(tused,string(round(Int,time()-startt)))    		
		push!(merrorlsls,merrorls)
    end	
	if isa(epsols_perind, AbstractVector) 
		offspringexcl = MagicReconstruct.del_error_indices(epsols_perind; tukeyfence,
			minoutlier=minoutlier,max_error=threshlikeparameters.peroffspringerror)						
		if !isempty(offspringexcl) 
			msg *= string(", #off_excl=",length(offspringexcl))
		end	
	end	
	errortuples = (epsf=epsfls,epso=epsols,  epso_perind = epsols_perind, seqerror = seqerrorls, allelebalancemean=allelebalancemeanls, 
		allelebalancedisperse=allelebalancedispersels, alleledropout = alleledropoutls)		
	if isdelmarker && (iteration > miditeration || !isimputefounder) && size(chrfhaplo,1)>2
		startt = time()		
		deltt = markerdelete_chr!(chrfhaplo,chroffgeno, popmakeup, priorprocess;
			errortuples..., offspringexcl, snporder,
			decodetempfile = imputetempfile,
			delsiglevel, israndallele, issnpGT, caldistance=isspacemarker,priorlength)
		ndel = length(deltt)
		ndel == 0 && (isdelmarker = false)
		msg *= string(", #del=",ndel)
		step_verbose && println(msg)
		push!(tused,string(round(Int,time()-startt)))
	end	
	if isordermarker && (iteration >= 2 || !isimputefounder) && size(chrfhaplo,1)>2		
		isneighbor = !isnothing(chrneighbor) && size(chrfhaplo,1)>10
		if isneighbor
			# update ordering, neighbor-based
			# logneighbor[i]: proptype, winsize, numprop, tused, apt, logl			
			logneighbor =updateorder_neighbor!(snporder,
				chrfhaplo,chroffgeno, popmakeup,priorprocess;
				chrneighbor, temperature, reversechr, orderactions = orderactions_neighbor, 
				errortuples...,  offspringexcl, maxwinsize=nbrmaxwin,
				israndallele,issnpGT,decodetempfile=imputetempfile)
			winsizels = [i[2] for i in logneighbor] 
			tusedls = [sum([i[4] for i in logneighbor])]						
			aptls = [i[5] for i in logneighbor]
			for i in eachindex(aptls)
				aptls[i] < minaccept && (nbrmaxwin[i] = max(3,round(Int,winsizels[i])))
			end
			aptls = [mean(aptls)]
			winsizels = [mean(winsizels)]
		else
			logneighbor  = []
			winsizels = []
			aptls = []
			tusedls = []
		end		
		logbook =updateorder!(snporder, chrfhaplo,chroffgeno, popmakeup,priorprocess; 
			errortuples..., offspringexcl, israndallele, issnpGT, slidewinsize, maxwinsize, orderactions, 
			decodetempfile=imputetempfile, reversechr, temperature)	
		push!(winsizels, mean([i[2] for i in logbook]))
		push!(tusedls, sum([i[4] for i in logbook]))
		push!(aptls, mean([i[5] for i in logbook]))
		isneighbor && (logbook = vcat(logneighbor,logbook))
		pushfirst!(logbook,["proptype", "win","nprop","tused","apt","logl"])
		logbook_order = permutedims(reduce(hcat,logbook))	
		winsizels2 = temperature <= 0.01 ? round.(winsizels,digits=1) : round.(Int, winsizels) 
		msg *= string(", win=",join(winsizels2, "_"))
		aptmsg = [isnan(i) ? NaN : round(Int,i) for i in 100*aptls]
		msg *= string(", apt=",join(aptmsg,"_"), "%")			
		push!(tused,join(round.(Int,tusedls),"_"))	
		step_verbose && println(msg)		
		bool = !any([isimputefounder,iscorrectfounder,isdelmarker])
		# && pdf(slidewinsize,2) >= 0.45   # pdf(truncated(Poisson(2),2,100), 2) ≈ 0.4556 		
		# chromosome is not revrse if isodd
		if bool && temperature == 0.0 && iteration >= miniteration_order && isodd(iteration)
			isordermarker = false
		end
	else
		logbook_order = nothing
	end
    if isspacemarker && (iteration >= 2 || !isimputefounder) 
        startt = time()
		markerspace_fun = spacebyviterbi ? markerspace_chr_viterbi! : markerspace_chr!				
		markerspace_fun(chrfhaplo,chroffgeno, popmakeup,priorprocess;
			errortuples...,  offspringexcl, snporder,reversechr,israndallele,issnpGT,
			decodetempfile = imputetempfile,temperature, itmax = iteration <= 4 ? 20 : 10)
        bool = !any([isimputefounder,iscorrectfounder,isdelmarker])
		# if bool
		if iteration > miditeration && temperature <= 0.25
            ndel3,msg_trim = trimchrend!(priorprocess; trimcm,trimfraction)			
			ndel3 > 0 && (msg *= string(", #del_trim=",ndel3, "(", msg_trim,")"))
        else
            ndel3 = 0
        end
        pri1=first(values(priorprocess))
        chrlen = 100*sum(pri1.markerdeltd[pri1.markerincl][1:end-1])
        bool = !any([isimputefounder,isdelmarker,iscorrectfounder,isinfererror,isordermarker])		
		if bool && temperature==0.0 && ndel3 == 0 && length(chrlenls) >=2
			diff_chrlen = min(abs.(chrlen .- chrlenls[end-1:end])...)
			nactionls = reverse(sum.(actionlsls))
			nstuck = findfirst(x->x>1, nactionls) -1 
            if diff_chrlen<=1 || diff_chrlen/chrlen < 0.005 || (nstuck >=5 && chrlen < chrlenls[end])
                isspacemarker = false
            end
        end
        push!(chrlenls,chrlen)
        msg *= string(", l=",round(Int,chrlen),"cM")
        step_verbose && println(msg)
        push!(tused,string(round(Int,time()-startt)))
    end	
	push!(actionlsls, [isimputefounder, iscorrectfounder, isinfererror, isinferseq, isdelmarker, isspacemarker,isordermarker])
	likeerrortuple = (epsfls=epsfls,epsols=epsols, epsols_perind=epsols_perind, seqerrorls=seqerrorls, 
		allelebalancemeanls = allelebalancemeanls, allelebalancedispersels = allelebalancedispersels, alleledropoutls = alleledropoutls)	
	likeerrortuple, offspringexcl, correctdf, tused, msg, logbook_order, upbyhalf
end


function save_impute_refine_chr!(magicgeno::MagicGeno, chr::Integer, chrfhaplo::AbstractMatrix,
	likeerrortuple::NamedTuple,offspringexcl::AbstractVector;
	isfounderinbred::Bool,isinfererror::Bool)
	(; epsfls, epsols, epsols_perind, seqerrorls,allelebalancemeanls,allelebalancedispersels, alleledropoutls) = likeerrortuple
	# set magicmped
	chrid = magicgeno.markermap[chr][1,:linkagegroup]		
	if !isnothing(epsols_perind)
		offinfo = magicgeno.magicped.offspringinfo
		colerr = string("peroffspringerror_",chrid)
		if colerr in names(offinfo)
			offinfo[!,colerr] .= epsols_perind
		else
			insertcols!(offinfo, size(offinfo,2)+1, colerr =>epsols_perind)
		end
		colexcl = string("offspringexcl_",chrid)
		offexcl = falses(size(offinfo,1))
		offexcl[offspringexcl] .= true
		if colexcl in names(offinfo)
			offinfo[!,colexcl] .= offexcl
		else
			insertcols!(offinfo, size(offinfo,2)+1, colexcl =>offexcl)
		end
	end
	# set chrfhaplo
	if isfounderinbred
		magicgeno.foundergeno[chr] .= Matrix{String}(chrfhaplo)
		magicgeno.markermap[chr][!,:founderformat] .=  "GT_haplo" 			
	else			
		magicgeno.foundergeno[chr] = Vector{String}.(permutedims(reduce(hcat,[Vector.(collect(partition(i,2))) for i in eachrow(chrfhaplo)])))
		magicgeno.markermap[chr][!,:founderformat] .= "GT_phased"			
	end
	# set error rates
	if isinfererror
		# keep errorrate as input if isinfererror = false
		magicgeno.markermap[chr][!,:foundererror] .= epsfls  
		magicgeno.markermap[chr][!,:offspringerror] .= epsols			
		magicgeno.markermap[chr][!,:seqerror] .= seqerrorls			
		magicgeno.markermap[chr][!,:allelebalancemean] .= allelebalancemeanls			
		magicgeno.markermap[chr][!,:allelebalancedisperse] .= allelebalancedispersels			
		magicgeno.markermap[chr][!,:alleledropout] .= alleledropoutls
	end				
	magicgeno
end

function get_winsizels(initwin::Real,niter::Integer)    
    x = initwin*1.0
    res = Vector{Float64}()
    push!(res,round(x,digits=1))  
	f(x) = x- sqrt(x)
    while true
        x = f(x)
        push!(res,max(2.0,round(x,digits=1)))
        x <= 2 && break  
    end
    if length(res) < niter    
        append!(res, 2*ones(niter-length(res)))
    end
    n = length(res)  
    winls = collect(range(initwin,2; length=n))
    map(max, res,winls)    
end

function get_likeerrortuple(magicgeno::MagicGeno,chr::Integer; 
    likeparameters::LikeParameters,
    isinfererror::Bool,
	repeatrun::Union{Nothing,Integer}=nothing, 
    io::Union{Nothing,IO}, 
    verbose::Bool)    
    liketargetls, epsf, epso, epso_perind, seqerror, allelebalancemean,allelebalancedisperse,alleledropout = MagicBase.extract_likeparameters(likeparameters)		
	chrmarkermap = magicgeno.markermap[chr]
    if isinfererror			
        epsfls = [ismissing(i) ? epsf : i for i in chrmarkermap[!,:foundererror]]
        epsols = [ismissing(i) ? epso : i for i in chrmarkermap[!,:offspringerror]]			
        seqerrorls = [ismissing(i) ? seqerror : i for i in chrmarkermap[!,:seqerror]]
        allelebalancemeanls = [ismissing(i) ? allelebalancemean : i for i in chrmarkermap[!,:allelebalancemean]]						
        allelebalancedispersels = [ismissing(i) ? allelebalancedisperse : i for i in chrmarkermap[!,:allelebalancedisperse]]			
        alleledropoutls = [ismissing(i) ? alleledropout : i for i in chrmarkermap[!,:alleledropout]]			        
        offinfo = magicgeno.magicped.offspringinfo
		chrid = chrmarkermap[1,:linkagegroup]
        colerr = string("peroffspringerror_",chrid)			
        if colerr in names(offinfo)				
            epsols_perind = [ismissing(i) ? epso_perind : i for i in offinfo[!,colerr]]
        else			
			noff = size(offinfo,1)
            epsols_perind = epso_perind*ones(noff)
        end
        merrls = [mean(epsfls),mean(epsols),mean(epsols_perind), mean(seqerrorls),mean(allelebalancemeanls),mean(allelebalancedispersels),mean(alleledropoutls)]			
		msg_run = isnothing(repeatrun) ? "" : string(", run=",repeatrun)
        msg = string("chr=", chrid, msg_run, ", initial <likeparameters> = ", round.(merrls,digits=4))
        printconsole(io,verbose, msg)
    else
        epsfls = epsf
        epsols = epso
        seqerrorls = seqerror
        allelebalancemeanls = allelebalancemean
        allelebalancedispersels = allelebalancedisperse
        alleledropoutls = alleledropout        
        epsols_perind = nothing
    end		
	likeerrortuple = (epsfls=epsfls,epsols=epsols, epsols_perind=epsols_perind, seqerrorls=seqerrorls, 
		allelebalancemeanls = allelebalancemeanls, allelebalancedispersels = allelebalancedispersels, alleledropoutls = alleledropoutls)
    liketargetls,likeerrortuple
end

function infer_juncdist(chrfhaplo::AbstractMatrix, chroffgeno::AbstractMatrix,
    popmakeup::AbstractDict,
    priorprocess::AbstractDict, priorspace::AbstractDict; 	
    epsf::Union{Real,AbstractVector},
    epso::Union{Real,AbstractVector},
	epso_perind::Union{Nothing, AbstractVector}, 
	seqerror::Union{Real,AbstractVector},
    allelebalancemean::Union{Real,AbstractVector},
    allelebalancedisperse::Union{Real,AbstractVector},
    alleledropout::Union{Real,AbstractVector},	
    snporder::AbstractVector,
    israndallele::Bool=true,     
	issnpGT::AbstractVector, 
    decodetempfile::AbstractString)
    # calculate chr_viterbi
    MagicReconstruct.hmmdecode_chr(chrfhaplo,chroffgeno,popmakeup,priorprocess;
        epsf,epso, epso_perind, seqerror,allelebalancemean,allelebalancedisperse,alleledropout,
		hmmalg="viterbi", decodetempfile, israndallele,issnpGT,snporder)
    nstate, nfgl = MagicReconstruct.hmm_nstate_nfgl(popmakeup)	
    if nstate == nfgl
        ancestralstate = [[i,i] for i in priorspace["haploindex"]]
    else
        nstate == nfgl^2 || @error "inconsistent number of states"
        ancestralstate = MagicBase.prior_diploindex(nfgl)
    end 
    resviterbi = last.(last(MagicReconstruct.read_decodefile(decodetempfile)))
    #  calculate juncdist
    avgnrecom = mean([length(MagicBase.splitindex(i[i.>0]))-1 for i in resviterbi])    
	pri1=first(values(priorprocess))
    chrlen = sum(pri1.markerdeltd[pri1.markerincl][1:end-1])
    mapexpansion = round(avgnrecom/chrlen,digits=4)
    nfounder = length(unique(reduce(vcat, [i["founder"] for i in values(popmakeup)])))
    juncdist = JuncDist(;nfounder,mapexpansion)
    juncdist
end

function up_fhaplopp!(fhaplosetpp::AbstractVector,correctdf::AbstractDataFrame)
	# println("corrrectdf=",correctdf)
	# println("fhaplosetpp[1][1:10]=",fhaplosetpp[1][1:10])
    for i in 1:size(correctdf,1)
        # [snp, p, oldfg,newfg,oldnerr,newnerr,nnonmiss]
        snp, p,oldg,newg = correctdf[i,1:4]
		msg = string("inconsisent founder correction. possible genotypes: ",
			fhaplosetpp[p][snp], ", oldg=",oldg, ",newg", newg, ",snp=", snp, ",p=",p)							
		if isa(oldg, AbstractVector) 
			if length(oldg) == 2
				if in("N",oldg)
					fhaplosetpp[p][snp] = [["1","1"],["1","2"],["2","1"],["2","2"],["1","N"],["2","N"],["N","1"],["N","2"],["N","N"]]
				else
					fhaplosetpp[p][snp] = [["1","1"],["1","2"],["2","1"],["2","2"]]
				end
			else
				@warn msg
			end
		else
			if oldg == "N" 
				fhaplosetpp[p][snp] = ["1","2","N"]
			else
				fhaplosetpp[p][snp] = ["1","2"]
			end
		end
    end
end

function index_founder2progeny(magicped::MagicPed)
    nfounder = size(magicped.founderinfo,1)
    frule = Dict(magicped.founderinfo[!,:individual] .=> 1:nfounder)
    offrule = Dict(magicped.offspringinfo[!,:individual] .=> 1:size(magicped.offspringinfo,1))
    founder2off = MagicBase.get_founder2offspring(magicped)
    res = Vector{Vector{Int}}(undef, nfounder)
    for (f,offs) in founder2off
        f2 = get(frule,f,nothing)
        res[f2] = [get(offrule,i,nothing) for i in offs]
    end
    res
end

function check_issnpGT(chroffgeno::AbstractMatrix,
    issnpGT::AbstractVector,errlab::AbstractString)
    b = isa.(chroffgeno[.!issnpGT,1],AbstractVector)
    if !all(b)
        @error string(errlab*": inconsistent issnpGT, chroffgeno[.!issnpGT,1]=",chroffgeno[.!issnpGT,1])
    end
end