

function get_subld(dosegeno::AbstractMatrix,snprange::AbstractRange,
    minlodsave::Real, minldsave::Real,
    show_progress::Bool, outfile::AbstractString,verbose::Bool)
    startt = time()
    nmarker = size(dosegeno,2)
    open(outfile,"w") do io
        progress = Progress(length(snprange);
            enabled = show_progress,
            desc=string("worker=",myid()," LD analyzing..."))
        for snp in snprange                
            nonmiss = view(dosegeno,:,snp) .> 0
            any(nonmiss) || continue
            geno_snp = dosegeno[nonmiss,snp]            
            length(unique(geno_snp)) <=1 && continue
            for snp2 in snp+1:nmarker
                # calculate LD = r^2            
                geno_snp22 = dosegeno[nonmiss, snp2]
                nonmiss2 = geno_snp22 .> 0      
                any(nonmiss2) || continue    
                geno_snp2 = geno_snp22[nonmiss2]
                length(unique(geno_snp2)) <=1 && continue
                geno_snp1 = geno_snp[nonmiss2]
                r = cor(geno_snp1, geno_snp2)
                isnan(r) && continue
                r2 = r^2
                r2 <= minldsave && continue
                # calculate independence LOD score
                crosstab = cal_crosstab(geno_snp1,geno_snp2)
                isnothing(crosstab) && continue
                gstat,df = cal_gstat_df(crosstab)
                df ≈ 0.0 && continue
                lod = cal_indeplod(gstat,df)
                lod <= minlodsave && continue
                # save results
                res = (snp,snp2,round(r2,digits=4),round(lod,digits=4))
                write(io, join(res,","),"\n")
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

function get_subld(dosegeno::AbstractSparseMatrix,snprange::AbstractRange,
    minlodsave::Real, minldsave::Real,
    show_progress::Bool, outfile::AbstractString,verbose::Bool)
    startt = time()
    nmarker = size(dosegeno,2)
    open(outfile,"w") do io
        progress = Progress(length(snprange);
            enabled = show_progress,
            desc=string("worker=",myid()," LD analyzing..."))
        for snp in snprange    
            nonmiss,geno_snp = findnz(dosegeno[:,snp])
            isempty(nonmiss) && continue
            length(unique(geno_snp)) <=1 && continue
            for snp2 in snp+1:nmarker
                # calculate LD = r^2            
                nonmiss2,geno_snp2 = findnz(dosegeno[nonmiss,snp2])
                isempty(nonmiss2) && continue    
                length(unique(geno_snp2)) <=1 && continue
                geno_snp1 = geno_snp[nonmiss2]
                r = cor(geno_snp1, geno_snp2)
                isnan(r) && continue
                r2 = r^2
                r2 <= minldsave && continue
                # calculate independence LOD score
                crosstab = cal_crosstab(geno_snp1,geno_snp2)
                isnothing(crosstab) && continue
                gstat,df = cal_gstat_df(crosstab)
                df ≈ 0.0 && continue
                lod = cal_indeplod(gstat,df)
                lod <= minlodsave && continue
                # save results
                res = (snp,snp2,round(r2,digits=4),round(lod,digits=4))
                write(io, join(res,","),"\n")
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

function cal_crosstab(geno_snp1::AbstractVector,geno_snp2::AbstractVector; 
    doseerror::Real=0.005)    
    #  analyze snp1
    rowname = unique(geno_snp1)
    length(rowname)<=1 && return nothing
    rowpos = [geno_snp1 .== i for i in rowname]
    nerror = min(max(1,length(geno_snp1)*doseerror),5) 
    b = sum.(rowpos) .> nerror # remove genotypes with frequency <= nerror
    rowname = rowname[b]
    rowpos = rowpos[b]
    length(rowname)<=1 && return nothing
    #  analyze snp2
    colname = unique(geno_snp2)
    length(colname)<=1 && return nothing
    cc = [[sum(geno_snp2[pos] .== j) for j = colname] for pos = rowpos]
    b = sum(cc) .> nerror
    cc2 = hcat(cc...)'
    crosstab = cc2[:,b]
    crosstab
end

function cal_gstat_df(counts::AbstractMatrix)
    rowsum = sum(counts, dims=2)
    colsum = sum(counts, dims=1)
    isempty(counts) && return zeros(2)
    ncol = sum(colsum .> 0)
    nrow = sum(rowsum .> 0)
    if ncol<=1 || nrow<=1
        zeros(2)
    else
        nn = sum(colsum)
        expect = kron(rowsum,colsum)/nn
        nonzero  = counts .!= 0
        ls =  log.(counts[nonzero] ./ expect[nonzero])
        gstat = 2*sum(counts[nonzero] .* ls)
        df  = (ncol -1) * (nrow -1)
        [gstat,df]
    end    
end

function cal_indeplod(gstat::Real,df::Real)
    # df ≈ 0.0 && return 0.0
    lcdf = logcdf(Chisq(df),gstat)
    lod = invlogcdf(Chisq(1),lcdf)
    # when g is too large, lod could be NaN or Inf
    (isinf(lod) || isnan(lod)) && (lod = gstat)
    lod *= log10(ℯ)/2
    lod
end
   