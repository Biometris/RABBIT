
function binning(genofile::AbstractString, 
    pedinfo::Union{MagicBase.JuncDist,AbstractString};
    isdepmodel::Bool=false,    
    formatpriority::AbstractVector=["GT","AD"],
    threshcall::Real = isdepmodel ? 0.95 : 0.9,   
    isfounderinbred::Bool=true,
    baseerror::Real=0.001,
    markerthin::Integer=1,
    binshare::Real=0.5,
    commentstring::AbstractString="##",
    isparallel::Bool=true,
    workdir::AbstractString=pwd(),
    outstem::AbstractString="outstem",
    logfile::Union{AbstractString,IO}= string(outstem,"_binning.log"),
    verbose::Bool=true)
    starttime = time()
    logio = MagicBase.set_logfile_begin(logfile, workdir, "binning"; verbose)
    msg = string("list of args: \n",
        "genofile = ", genofile, "\n",
        "pedinfo = ", pedinfo, "\n",
        "isdepmodel = ", isdepmodel, "\n",
        "formatpriority = ", formatpriority, "\n",
        "threshcall = ", threshcall, "\n",
        "isfounderinbred = ", isfounderinbred, "\n",
        "baseerror = ", baseerror, "\n",
        "markerthin = ", markerthin, "\n",
        "binshare = ", binshare, "\n",
        "commentstring = ", commentstring, "\n",
        "isparallel = ", isparallel, isparallel ? string("(nworkers=",nworkers(),")") : "", "\n",
        "workdir = ",workdir,"\n",
        "outstem = ", outstem, "\n",
        "logfile = ",logfile,"\n",
        "verbose = ",verbose)
    printconsole(logio,verbose,msg)
    ispath(workdir) || @error string(workdir, " is not a directory")
    magicgeno = formmagicgeno(genofile, pedinfo;
        formatpriority,isfounderinbred, commentstring,workdir)
    merge_chromosome!(magicgeno)
    nsnp = size(only(magicgeno.markermap),1)
    MagicBase.submagicgeno!(magicgeno, snpsubset=1:markerthin:nsnp)
    MagicBase.info_magicgeno(magicgeno;io=logio,verbose)
    MagicBase.rawgenoprob!(magicgeno; baseerror,isfounderinbred,
        isoffspringinbred = isdepmodel)
    MagicBase.rawgenocall!(magicgeno; isfounderinbred, callthreshold = threshcall)        
    # set A2 as minor aallele
    tused = @elapsed consistent_allele!(magicgeno)
    printconsole(logio,verbose,string("tused=",round(tused,digits=1),"s in setting consistent minor allele"))
    tused = @elapsed ishalfmiss, codegeno = get_codegeno(magicgeno;isdepmodel);
    printconsole(logio,verbose,string("tused=",round(tused,digits=1),"s in coding genotypes"))
    markermap = deepcopy(only(magicgeno.markermap))
    snpidls = markermap[:,:marker]
    nsnp = length(snpidls)
    nwork = nworkers()
    mem1 = round(Int, memoryuse()/10^6)
    magicgeno = nothing
    GC.gc()
    mem2 = round(Int, memoryuse()/10^6)
    msg = string("memoryuse=",mem1,"|", mem2,"MB")  
    printconsole(logio,verbose,msg)  
    if isparallel && nwork>1
        if nsnp <= 10*nwork
            nseg = 1
        else
            nseg0 = round(Int, nsnp/(1000*nwork))
            nseg = (nseg0+1)*nwork
        end
    else
        nseg = nsnp <= 1000 ? 1 : div(nsnp,1000)
    end
    # snprgls = MagicLD.getsnprg(nsnp,nseg)    
    snprgls = collect(partition(1:nsnp,ceil(Int,nsnp/nseg)))    
    # startadj = time()
    msg = string(length(snprgls), " sets of snpranges")
    if nwork > 1
        msg *= string(", distributed over ", nwork, " workers")        
    end
    printconsole(logio,verbose,msg)
    filels = [joinpath(workdir,string(outstem,"_binning_temporary",i,".txt")) for i=1:length(snprgls)]
    startparallel = time()
    if isparallel && nwork>1
        msgls = pmap((x,y)->pairwisedupe(codegeno, ishalfmiss, x, y,binshare,false,verbose), snprgls,filels)
    else
        msgls = map((x,y)->pairwisedupe(codegeno, ishalfmiss, x, y,binshare,true,verbose), snprgls,filels)
    end
    mem1 = round(Int, memoryuse()/10^6)
    msg = string("tused = ", round(time()-startparallel,digits=1),"s for all markers, mem=",mem1,"MB")
    printconsole(logio, verbose,msg)
    for msg = msgls
        write(logio,string(msg,"\n"))
    end
    flush(logio)
    adjfile = string(outstem,"_binning_adjacency.csv.gz")
    adjfile2 = getabsfile(workdir,adjfile)
    GZip.open(adjfile2,"w") do outio        
        initial = "RABBIT"
        delim = ','
        nmissingls = count(==(0), codegeno, dims=1)[1,:]
        markerinfo = DataFrame(markerno = 1:length(snpidls), marker=snpidls,nmissing = nmissingls)
        MagicBase.appenddf(outio, markerinfo; delim,initial,dfname="markerinfo")
        write(outio,initial,delim,"adjacency\n")
        write(outio, "marker1, marker2\n")
        for file = filels
            for line in eachline(file,keep=true)
                 write(outio,line)
             end
             rm(file)
        end
    end
    # msg = string("calculating adjacency, tused = ", round(time()-startadj,digits=1),"s")
    msg = "Start binning from adjacency; it might take a while..."
    printconsole(logio,verbose,msg)
    tused = @elapsed binls, snpmiss = binning_dupe(adjfile2; delim=',', workdir)
    msg = string("binning from adjacency, tused = ", round(tused,digits=1),"s")
    printconsole(logio,verbose,msg)
    # codegeno rule: Dict(["1"=>1,"11"=>1, "2"=>2, "22"=>2,"12"=>3,"1N"=>4,"2N"=>5])  with default "N"=>0 and "NN"=>0
    # change into: Dict(["N"=>0.0,"NN"=>0.0,"1N"=>0.5,"2N"=>0.5]) and others("1","2","12","11","22") to 1.0    
    # codegeno[ismissing.(codegeno)] .= 0.0
    codegeno[codegeno .>= 4] .= 0.5 # half observed 
    codegeno[codegeno .>= 1] .= 1 # full osbvered    
    bindf = binning_share(codegeno, markermap, binls, snpmiss; binshare)
    binfile = string(outstem,"_binning.csv.gz")
    binfile2 = joinpath(workdir,binfile)
    GZip.open(binfile2,"w") do io
        descripls = ["marker, marker ID", 
            "linkagegroup, linkage group ID", 
            "poscm, marker position in centiMorgan",
            "physchrom, physical map chromosome ID for the marker",
            "physposbp, physical map position (in base pair) for the marker",
            "binno, marker bin index",
            "represent, same as binno if the marker is representative and 0 otherwise"
        ]
        for i in eachindex(descripls)
            write(io, string("##col_",i, ", ", descripls[i],"\n"))
        end
        CSV.write(io, bindf; delim=',', header=true,missingstring="NA", append=true)
    end    
    ubinno = unique(bindf[!,:binno])
    ndel = sum(ubinno .< 0)
    msg = string("#markers = ",size(snpmiss,1),
        ", #del_markers = ",  ndel,
        ", #bins = ", length(ubinno)-ndel)
    printconsole(logio,verbose,msg)
    msg = string("adjacency file: ", adjfile)
    printconsole(logio,verbose,msg)
    msg = string("binning file: ", binfile)
    printconsole(logio,verbose,msg)
    tused = round(time()-starttime,digits=1)
    MagicBase.set_logfile_end(logfile, logio, tused,"binning"; verbose)
    binfile2, adjfile2
end

function get_codegeno(magicgeno::MagicGeno;isdepmodel::Bool=false)    
    snpgeno = permutedims(hcat(only(magicgeno.foundergeno),only(magicgeno.offspringgeno)))
    gset = setdiff(unique(snpgeno),["N","NN"])
    # rule = Dict(["1"=>1,"11"=>1, "2"=>2, "22"=>2,"12"=>3,"1N"=>4,"2N"=>5])
    codetype = Float16
    rule = Dict(["1"=>codetype(1),"11"=>codetype(1), "2"=>codetype(2), "22"=>codetype(2),
        "12"=>codetype(3),"1N"=>codetype(4),"2N"=>codetype(5)])
    d = setdiff(gset,keys(rule))
    isempty(d) || @error string("unexpected genotypes: ",d)         
    if isdepmodel
        delete!(rule,"12")
        delete!(rule,"21")
        setdiff!(gset, ["12","21"])    
    end
    ishalfmiss = !isempty(intersect(gset,["1N","2N"]))
    if !ishalfmiss 
        delete!(rule,"1N")
        delete!(rule,"2N")
    end    
    codemiss = codetype(0)    
    res = [get(rule,i,codemiss) for i in snpgeno]
    ishalfmiss, res
end

function readdupefile(dupefile::AbstractString;
    delim::AbstractChar=',',
    commentstring::AbstractString = "##", 
    workdir::AbstractString=pwd())
    dupefile2 = getabsfile(workdir,dupefile)    
    res= MagicBase.readmultitable(dupefile2;
        delim,
        missingstring=["missing"], # keep NA 
        commentstring
    )
    markerinfo = res["markerinfo"]
    markers = string.(markerinfo[!,:marker])    
    nmissingls = markerinfo[!,:nmissing]
    df = res["adjacency"]
    nsnp = length(markers)
    iils, jjls = df[!,1], df[!,2]
    iils, jjls = vcat(iils, jjls, 1:nsnp), vcat(jjls, iils, 1:nsnp)
    vvls = ones(length(iils))
    adjdupe  = sparse(iils,jjls, vvls, nsnp,nsnp) # default combine option is +
    inputno = 1:length(markers)
    snpmiss = DataFrame("marker"=>markers,"inputno"=>inputno, "nmissing"=>nmissingls)
    adjdupe, snpmiss
end

function binning_dupe(dupefile::AbstractString;    
    delim::AbstractChar=',',
    workdir::AbstractString=pwd())
    adjdupe, snpmiss = readdupefile(dupefile; delim, workdir)
    binls = binning_dupe(adjdupe,snpmiss[!,:nmissing])    
    binls, snpmiss
end

function binning_dupe(adjdupe::AbstractMatrix,nmissingls::AbstractVector)
    nsnp = length(nmissingls) 
    size(adjdupe,1) == nsnp || @error "inconsistent #markers"
    gdupe = SimpleGraph(adjdupe)
    vvls = partition_graph(gdupe)    
    res = ThreadsX.map(vvls) do x
        MagicMap.binning_dupe_sub(gdupe,nmissingls,x)
    end
    binls = reduce(vcat, res)
    binls = binls[sortperm(mean.(binls))]
    @assert sort(reduce(vcat,binls)) == 1:nsnp # all markers are binned
    binls
end

function binning_dupe_sub(gdupe::SimpleGraph, nmissingls::AbstractVector, vv::AbstractVector)
    g, vmap = induced_subgraph(gdupe,vv)        
    # @assert vmap == vv
    weight = -1 .* view(nmissingls, vmap)
    ls = partition_graph(g, weight)    
    [vv[i] for i in ls]
end

function binning_share(codegeno::AbstractMatrix, markermap::DataFrame,
    binls::AbstractVector,
    snpmiss::AbstractDataFrame; binshare::Real=0.5)
    markermap[!,:marker] == snpmiss[!,:marker] || @error "inconsistent markers"
    chrls = markermap[!,:linkagegroup] 
    poscmls = markermap[!,:poscm] 
    physchromls = markermap[!,:physchrom]         
    physposbpls = markermap[!,:physposbp]         
    nsnp = size(snpmiss,1)
    binno = zeros(Int,nsnp)
    represent = zeros(Int, nsnp)
    unbinned = Vector{Int}()
    for i in eachindex(binls)
        bin = binls[i]
        if length(bin)==1
            represent[bin] .= i
            binno[bin] .= i
        else
            snp = bin[argmin(snpmiss[bin,:nmissing])]
            represent[snp] = i
            restsnp = setdiff(bin, [snp])
            pos1 = codegeno[:,snp] .== 1.0
            share1 = sum(codegeno[pos1,restsnp], dims=1)[1,:] 
            pos05 = codegeno[:,snp] .== 0.5
            share05 = 0.5*sum(sign.(codegeno[pos05,restsnp]), dims=1)[1,:] 
            nshare = share1 .+ share05
            fshare = nshare ./ (sum(pos1) + 0.5*sum(pos05))            
            b = @. fshare > binshare && nshare >= 20
            keepsnp = restsnp[b]
            binno[keepsnp] .= i
            binno[snp] = i
            push!(unbinned,setdiff(restsnp,keepsnp)...)
        end
    end
    if !isempty(unbinned)
        # negative ids  denote marker deletion
        id = -1 * (1:length(unbinned))
        binno[unbinned] .= id
        represent[unbinned] .= id
    end
    bindf = DataFrame("marker"=>snpmiss[!,:marker],
        "linkagegroup"=>chrls,
        "poscm"=>poscmls,
        "physchrom" => physchromls,
        "physposbp"=>physposbpls,
        "binno"=>binno,
        "represent"=>represent)
    bindf
end

# maxdegree::Integer = 20, maxit::Integer = 200
function findcliques(g::AbstractGraph,v::Integer;
    maxdegree::Integer = 20, maxit::Integer = 200)
    d = degree(g,v)
    d == 0 && return [[v]]
    if d <= maxdegree
        cliques = localcliques(g,v)
    else
        cliques = Vector{Vector{Int}}()
        maxlen = 0
        nstuck = 0
        for it  in 1:ceil(Int,maxit/10)
            nres = length(cliques)
            ls =[MagicMap.localclique(g,v) for i in 1:10]        
            lenls = length.(ls)
            b = lenls .>= maxlen
            append!(cliques,sort.(ls[b]))            
            unique!(cliques)    
            if nres == length(cliques) 
                nstuck += 1
            else
                nstuck = 0
                maxlen = max(maxlen, max(lenls[b]...))
            end        
            # println("it=",it, ",len=",length(cliques), ",nstuck=",nstuck)
            nstuck >= 10 && break
        end
        len = length.(cliques)
        maxlen = max(len...)
        cliques = cliques[len .== maxlen]
    end
    cliques
end

function localcliques(g::AbstractGraph,v::Integer)
    nbrls = [setdiff(neighbors(g,i),[i]) for i in 1:nv(g)]
    nbr = nbrls[v]
    isempty(nbr) && return [[v]]
    res = [(sort([v,i]),sort(intersect(nbr,nbrls[i]))) for i in nbr]
    while true
        all([isempty(i[2]) for i in res]) && break
        ls = [begin
            clique,  common = res[i]
            # intersect keeps order with the first arg
            [(sort(vcat(clique,[j])), intersect(common,nbrls[j])) for j in common]
        end for i in eachindex(res)]
        res = unique(reduce(vcat,ls))
        # println("nclique = ",length(res))
    end
    cliques = [i[1] for i in res]
    cliques
end

function localclique(g::AbstractGraph,v::Integer)
    nbrls = [setdiff(neighbors(g,i),[i]) for i in 1:nv(g)]
    clique = [v]
    common = nbrls[v]
    isempty(common) && return clique
    while true
        j = rand(common)
        clique = sort(vcat(clique,[j]))
        common = intersect(common,nbrls[j])
        isempty(common) && break
    end
    clique
end

function partition_graph(gdupe::Graph)
    cc = connected_components(gdupe)
    res = [if length(c) <= 2
        [c]
    else
        g, vmap = induced_subgraph(gdupe,c)
        lab = first(label_propagation(g))
        [vmap[lab .== i] for i  in unique(lab)]
    end for c in cc]
    reduce(vcat,res)
end

function partition_graph(g::AbstractGraph, weight::AbstractVector)
    res = Vector{Vector{Int}}()
    gg = g
    vvlab = 1:nv(g)
    while true
        v = argmax(weight[vvlab])
        cliques = findcliques(gg,v)
        clique = intersect(cliques...)        
        vlist = setdiff(1:nv(gg),clique)
        push!(res, vvlab[clique])
        isempty(vlist) && break
        gg = first(induced_subgraph(gg,vlist))
        vvlab = vvlab[vlist]
        # println("v=", v,", clique_size=",length(clique),",tused=",tused,"s")
    end
    # @assert sort(vcat(res...)) == 1:nv(g)
    res
end


# set minor allele = allele2, phasing negected
function consistent_allele!(magicgeno::MagicGeno)
    for chr in 1:length(magicgeno.markermap)
        snpmap = magicgeno.markermap[chr]
        d = setdiff(unique(snpmap[!,:founderformat]),["GT_haplo","GT_unphased"])
        if !isempty(d)
            @error string("TODO for founerformat = ",d)
            continue
        end
        d = setdiff(unique(snpmap[!,:offspringformat]),["GT_haplo","GT_unphased"])
        if !isempty(d)
            @error string("TODO for offspringformat = ",d)
            continue
        end
        geno = hcat(magicgeno.foundergeno[chr],magicgeno.offspringgeno[chr])
        # set unphased geno        
        geno[geno .== "21"] .= "12"
        geno[geno .== "N1"] .= "1N"
        geno[geno .== "N2"] .= "2N"
        # set a2 as minor allele
        b = [begin 
            str = join(i)
            n1 = count("1",str)
            n2 = count("2",str)
            n1 < n2
        end for i in eachrow(geno)]
        subgeno = view(geno,b,:)
        subgeno .= replace.(subgeno,"2"=>"1","1"=>"2")
        subgeno[subgeno .== "21"] .= "12"
        nf = size(magicgeno.foundergeno[chr],2)
        magicgeno.foundergeno[chr] .= geno[:,1:nf]
        magicgeno.offspringgeno[chr] .= geno[:,(nf+1):end]
    end
    magicgeno
end

# geno[m,i]: genotypes of individual i at marker m
# function calfreq_a2(geno::AbstractMatrix)
#     n1 = sum(geno .== "1",dims=2) + sum(geno .== "1N",dims=2)
#     n2 = sum(geno .== "2",dims=2) + sum(geno .== "2N",dims=2)
#     n11 = sum(geno .== "11",dims=2)
#     n12 = sum(geno .== "12",dims=2)
#     n22 = sum(geno .== "22",dims=2)
#     tot_a2 = 2*(n22 .+ n12 .+ n11) .+ n2 .+ n1
#     count_a2 = 2*n22 .+ n12 .+ n2
#     freq_a2 = count_a2 ./ tot_a2
#     freq_a2[:,1]
# end

function pairwisedupe(codegeno::AbstractMatrix, ishalfmiss::Bool, snprange::UnitRange,
    outfile::AbstractString,
    binshare::Real,
    show_progress::Bool,
    verbose::Bool)
    startt = time()
    nmarker = size(codegeno,2)
    open(outfile, "w") do io
        # desc=string("worker=",myid(), " binning...")        
        progress = Progress(length(snprange);        
            enabled = show_progress,
            desc=string("binning..."))        
        if ishalfmiss            
            for snp in snprange
                # coding: Dict(["N"=>0,"NN"=>0,"1"=>1,"11"=>1, "2"=>2, "22"=>2,"12"=>3,"1N"=>4,"2N"=>5])
                # [abs(i-j)*i*j for i in 1:5, j in 1:5]
                # compatible between "11" and "1N": |"11"-"1N"|*"11"*"1N" = 12
                # compatible between "12" and "1N": |"12"-"1N"|*"12"*"1N" = 12
                # compatible between "22" and "2N": |"22"-"2N"|*"22"*"2N" = 30
                # compatible between "12" and "2N": |"12"-"2N"|*"12"*"2N" = 30
                # no other cases resulting in 12 and 30                                
                inds = findall(view(codegeno,:,snp) .> 0)
                length(inds) < 50 && continue      
                count_nonmiss = [sum(sign.(col)) for col in eachcol(view(codegeno,inds, snp+1:nmarker))]                              
                snp2ls = snp .+ findall(count_nonmiss .>= binshare*length(inds))                                
                ndiff = sum(begin
                    ls = abs.(codegeno[i,snp] .- codegeno[i,snp2ls])
                    ls .*= (codegeno[i,snp] .* codegeno[i,snp2ls])
                    ls .*= (ls .!= 12) .* (ls .!= 30)
                    sign.(ls)
                end for i in inds)
                snp2ls = snp2ls[ndiff .== 0]
                if !isempty(snp2ls)                    
                    msg = join([string(snp,",",j)  for j in snp2ls],"\n")
                    write(io, msg, "\n")
                end
                next!(progress)
            end
        else
            for snp in snprange
                # coding: Dict(["N"=>0,"NN"=>0,"1"=>1,"11"=>1, "2"=>2, "22"=>2,"12"=>3,"1N"=>4,"2N"=>5])
                nonmiss = view(codegeno,:,snp) .> 0    
                sum(nonmiss) < 50 && continue    
                geno_snp = codegeno[nonmiss,snp]    
                count_nonmiss = [sum(sign.(col)) for col in eachcol(view(codegeno,nonmiss, snp+1:nmarker))]                              
                snp2ls = snp .+ findall(count_nonmiss .>= binshare*length(nonmiss))              
                # snp2ls = findall(sum(view(codegeno,nonmiss, :),dims=1)[1,:] .== sum(geno_snp))
                # setdiff!(snp2ls,[snp])
                for snp2 in snp2ls      
                    geno_snp22 = codegeno[nonmiss,snp2]
                    nonmiss2 = geno_snp22 .> 0
                    any(nonmiss2) || continue                
                    b = MagicMap.isvecequal(geno_snp[nonmiss2], geno_snp22[nonmiss2])        
                    if b
                        write(io, join([snp,snp2],","),"\n")
                    end  
                end
                flush(io)
                next!(progress)
            end
        end
        mem1 = round(Int, memoryuse()/10^6)
        GC.gc()
        mem2 = round(Int, memoryuse()/10^6)
        msg = string("worker=",myid(),", snps=",snprange, 
            ", t=", round(time()-startt,digits=1), "s",
            ", mem=",mem1,"|", mem2,"MB")
        verbose && (@info msg)  
        msg
    end
end


function isvecequal(v1::AbstractVector,v2::AbstractVector)
    for i in eachindex(v1)
        v1[i] == v2[i] || return false
    end
    true
end

