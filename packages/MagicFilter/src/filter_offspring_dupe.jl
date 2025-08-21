
function filter_offspring_dupe!(magicgeno::MagicGeno;
    model::AbstractString="jointmodel",
    isfounderinbred::Bool=true,
    threshcall::Real = 0.9,
    offspring_maxcorr::Real = 0.99,
    offspring_cutcorr::Real = 0.4,
	isparallel::Bool=false,
    outstem::AbstractString= "outstem",
    logfile::Union{AbstractString,IO} = outstem*"_purify_dupe.log",
    workdir::AbstractString=pwd(),
    verbose::Bool=true)
    starttime = time()
    offspring_maxcorr > 1.0 && return magicgeno
    logio = MagicBase.set_logfile_begin(logfile, workdir, "filter_offspring_dupe!"; verbose)
    model = MagicBase.reset_model(magicgeno.magicped,model;io=logio,verbose)
    msg = string("list of options: \n",
        "model = ", model, "\n",
        "isfounderinbred = ", isfounderinbred, "\n",
        "threshcall =", threshcall, "\n",
        "offspring_maxcorr = ", offspring_maxcorr, "\n",
		"offspring_cutcorr = ", offspring_cutcorr, "\n",
		"isparallel = ", isparallel, isparallel ? string("(nworkers=",nworkers(),")") : "", "\n",
        "workdir = ",workdir,"\n",
        "outstem = ", outstem,"\n",
        "logfile = ",logfile,"\n",
        "verbose = ",verbose)
    MagicBase.printconsole(logio,verbose,msg)
    MagicBase.setunphasedgeno!(magicgeno)
    MagicBase.info_magicgeno(magicgeno;io=logio,verbose)
	# TOCheck: markermap is required in merge_missprogeny!
    missgeno =  MagicBase.split_missprogeny!(magicgeno)
    nsnpmiss = sum(size.(missgeno.markermap,1))
    if nsnpmiss > 0
        msg = string("remove (and insert afterwards) ", nsnpmiss, " markers with all offspring genotypes being missing")
        printconsole(logio,verbose,msg)
    end
    # find offspring duplicate
	tused = @elapsed  begin 
        dosegeno = MagicBase.get_dosegeno(magicgeno; baseerror=0.001, callthreshold=threshcall, iscalling=true, isdepmodel = model=="depmodel")                    
        mem1 = round(Int, memoryuse()/10^6)
        GC.gc()
        mem2 = round(Int, memoryuse()/10^6)        
    end
    msg = string("get_dosegeno, tused=", round(tused,digits=1), "s, mem=",mem1,"|",mem2,"MB")
    printconsole(logio, verbose, msg)
	nwork = nworkers()
	noffspring = size(dosegeno,2)
    if noffspring < 100 
		nseg = 1 
	else
		nseg =  (1 + div(noffspring,1000*(1+div(nwork,10)))) * nwork		
	end
    offrng = get_offspring_rg(noffspring,nseg)    
    msg = string(length(offrng), " sets of offspring ranges")
    nwork > 1 && (msg *= string(", distributed over ", nwork, " workers"))
    printconsole(logio,verbose,msg)   
    filels = [joinpath(workdir,string(outstem*"_filter_offspring_dupe_temporary",i,".txt")) for i=1:length(offrng)]
    startt = time()   
	if isparallel && length(offrng)>1
        show_progress = false
        msgls = pmap((x,y)->get_offspring_corr(dosegeno, x, offspring_cutcorr, show_progress,y,verbose), offrng,filels)
    else
        show_progress = true
        msgls = map((x,y)->get_offspring_corr(dosegeno, x, offspring_cutcorr, show_progress,y,verbose), offrng,filels)
    end
    for msg = msgls
        write(logio,string(msg,"\n"))
    end
    flush(logio)
    mem = round(Int, memoryuse()/10^6)
    msg = string("tused=",round(time()-startt, digits=1), "s for all offspring, mem=",mem, "MB")
    printconsole(logio, verbose, msg)
    outcorrfile = string(outstem,"_ind_correlation.csv.gz")    
    outcorrfile2 = getabsfile(workdir,outcorrfile)
    GZip.open(outcorrfile2,"w") do outio        
		write(outio, "RABBIT,offspringinfo\n")       
		CSV.write(outio, magicgeno.magicped.offspringinfo; header=true,delim=",",append=true)
		write(outio, "RABBIT,pairwisecorr\n")       
        write(outio, "offspring1,offspring2,correlation\n")       
        for file in filels
            if isfile(file)
                write(outio,read(file))
                flush(outio)
            end
        end        
    end	
    msg = string("pairwise offspring correlations: ", outcorrfile)
    printconsole(logio,verbose,msg)
    rm.(filels;force=true)
	relate, offinfo = read_offspring_corr(outcorrfile; workdir)
	offidls = offinfo[!,:individual]
    if nnz(relate) == 0
		msg = "no offspring deleted"
		printconsole(logio,verbose,msg)
		if nsnpmiss > 0
			MagicBase.merge_missprogeny!(magicgeno,missgeno; isspacemarker = false)
		end
	else		
        if length(offidls) < 2*10^4            
			relate2 = Matrix(relate)			
            gheat = heatmap(relate2,
                title = string("reset corr = 0 if |corr| < ",offspring_cutcorr),
                xlabel="Offspring index",
                ylabel="Offspring index",                
				color = cgrad(:bluesreds),
				size = (900,900), 
				left_margin=20Plots.mm,
				right_margin=10Plots.mm,				
        		bottom_margin=10Plots.mm,        
            )
            outfile = outstem*"_ind_correlation.png"
            savefig(gheat,joinpath(workdir,outfile))
        end
        # delete duplicates
		if offspring_maxcorr < 1.0
			nnonmiss = sum(dosegeno .> 0,dims=1)[1,:]			
			Is, Js, Vs = findnz(relate)
			b = Vs .> offspring_maxcorr
			off1 = Is[b]
			off2 = Js[b]			
			offdells = []
			for i in eachindex(off1,off2)
                off1[i] == off2[i] && continue
				w1, w2 = nnonmiss[off1[i]], nnonmiss[off2[i]]
				if w1 < w2
					push!(offdells,off1[i])
				elseif w1 > w2
					push!(offdells,off2[i])
				else
					off1[i] < off2[i] ? push!(offdells,off1[i]) : push!(offdells,off2[i])
				end
			end
			union!(offdells)
	        offkeep = collect(1:length(offidls))
	        if !isempty(offdells) 
                dells = reduce(vcat,offdells)
                # b = [occursin(r"_virtualoffspring$", offidls[i]) for i in dells]
                # deleteat!(dells,b)
                setdiff!(offkeep,dells)
            end
		end
		if nsnpmiss > 0
			MagicBase.merge_missprogeny!(magicgeno,missgeno; isspacemarker=false)
		end
        if offspring_maxcorr < 1.0            
			del_offspring!(magicgeno,offkeep;io=logio,verbose)
		end
    end
    # MagicBase.info_missing(magicgeno;io=logio,verbose)
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, logio, tused,"purify_dupe"; verbose)
    magicgeno
end

function get_offspring_rg(noffspring::Integer,nseg::Integer)
    ls = accumulate(+,noffspring-1:-1:1)
    ls2 = round.(Int,(1:(nseg-1)) .* (last(ls)/nseg))
    ls3 = [findfirst(x->x>i, ls) for i=ls2]
    pushfirst!(ls3,1)
    push!(ls3,noffspring)
    ls3 = unique(ls3)
    res = [ls3[i]:(ls3[i+1]-1) for i=1:length(ls3)-1]
    reverse(res)
end

function get_offspring_corr(dosegeno::AbstractMatrix,offrange::AbstractRange,
    offspring_cutcorr::Real, 		
    show_progress::Bool, outfile::AbstractString,verbose::Bool)
    startt = time()
	noffspring = size(dosegeno,2)
	open(outfile,"w") do io
        progress = Progress(length(offrange);
            enabled = show_progress,
            desc=string("worker=",myid()," Offspring correlation..."))
        for off in offrange                
			# missing is denoted by 0.0
            nonmiss = view(dosegeno,:,off) .> 0 
            geno_off = dosegeno[nonmiss,off]
            any(nonmiss) || continue
            length(unique(geno_off)) <=1 && continue			
            for off2 in off+1:noffspring
                # calculate LD = r^2            
                geno_off22 = dosegeno[nonmiss, off2]
                nonmiss2 = geno_off22 .> 0      
                any(nonmiss2) || continue    
                geno_off2 = geno_off22[nonmiss2]
                length(unique(geno_off2)) <=1 && continue
                geno_off1 = geno_off[nonmiss2]
                r = cor(geno_off1, geno_off2)
                isnan(r) && continue                
                r < offspring_cutcorr && continue                
                # save results
                res = (off,off2,round(r,digits=4))
                write(io, join(res,","),"\n")
            end
            flush(io)
            next!(progress)
        end
    end    
    mem1 = round(Int, memoryuse()/10^6)
    GC.gc()
    mem2 = round(Int, memoryuse()/10^6)
    msg = string("worker=",myid(),", offspring=",offrange, 
        ", t=", round(time()-startt,digits=1), "s",
        ", mem=",mem1,"|", mem2,"MB")
    verbose && (@info msg)    
    msg	
end


function read_offspring_corr(indcorrfile::AbstractString;
    delim::AbstractChar=',',
    workdir::AbstractString=pwd())
    indcorrfile2 = getabsfile(workdir, indcorrfile)
    resdict= MagicBase.readmultitable(indcorrfile2; delim)
    offinfo = resdict["offspringinfo"]
    noff = size(offinfo,1)
    corr = resdict["pairwisecorr"]
    relate = sparse(corr[!,1],corr[!,2],corr[!,3],noff,noff)
    relate .+= relate'
    relate, offinfo
end




function cal_off_relate(calledgeno::AbstractMatrix; offspring_cutcorr=0.4)
    gset = unique(calledgeno)
	rule = Dict("11"=>0.0,"12"=>1.0,"21"=>1.0,"22"=>2.0,"NN"=>missing,
	    "1N"=>0.5,"N1"=>0.5,
		"2N"=>1.5,"N2"=>1.5)
	# wrong correlation if using Float16 in case of large dataset	
	codegeno = Matrix{Union{Missing,Float64}}(undef,size(calledgeno)...)
	for g in gset
	    codegeno[calledgeno .== g] .= rule[g]
	end
	noff = size(codegeno,2)
	relate = spzeros(noff,noff)	
	nthread = Threads.nthreads()    
	progress = Progress(noff;desc="calculating offspring correlation...")       
	for i in 1:noff
		relate[i,i] = 0.5
		bi = .!ismissing.(view(codegeno,:,i) )
		if !any(bi) 
			next!(progress)
			continue
		end
		geno_i = codegeno[bi,i]		
		ls = zeros(noff-i)				
		basesize = max(1,div(length(ls),3*nthread))
		ThreadsX.foreach(eachindex(ls);basesize) do j
			geno_j = codegeno[bi,i+j]
			b = .!ismissing.(geno_j)
			if any(b) 
				r = abs(cor(geno_i[b],geno_j[b]))
				r >= offspring_cutcorr && (ls[j] = r)			
			end
		end
		relate[i,(i+1):noff] .= ls
		next!(progress)
	end
    nnonmiss = [begin
        c = 1.0*sum(.!ismissing.(i))
        c -= sum(skipmissing(i) .== 0.5)/2
        c -= sum(skipmissing(i) .== 1.5)/2
        c
    end for i in eachcol(codegeno)]
    relate, nnonmiss
end

function cal_off_relate2(calledgeno::AbstractMatrix; offspring_cutcorr=0.4)
	gset = unique(calledgeno)
	rule = Dict("11"=>0.0,"12"=>1.0,"21"=>1.0,"22"=>2.0,
	    "NN"=>missing,
	    "1N"=>missing,"N1"=>missing,"2N"=>missing,"N2"=>missing)
	codegeno = Matrix{Union{Missing,Float64}}(undef,size(calledgeno)...)
	for g in gset
	    codegeno[calledgeno .== g] .= rule[g]
	end
	noff = size(codegeno,2)
	relate = spzeros(Float64,noff,noff)
	@showprogress 1 "calculating offspring_corr..." for i in 1:noff
	    gi = view(codegeno,:,i)
	    bi = .!ismissing.(gi)
	    gi = gi[bi]
	    subgeno = view(codegeno, bi, i+1:noff)
	    res = [begin
	        d = col .- gi
	        mean(d[.!ismissing.(d)] .== 0.0)
	    end for col in eachcol(subgeno)]
	    res[res .< offspring_cutcorr] .= 0.0
	    relate[i,i+1:noff] .= res
	end
    nnonmiss = [begin
        c = 1.0*sum(.!ismissing.(i))
        # c -= sum(skipmissing(i) .== 0.5)/2
        # c -= sum(skipmissing(i) .== 1.5)/2
        c
    end for i in eachcol(codegeno)]
    relate, nnonmiss
end
