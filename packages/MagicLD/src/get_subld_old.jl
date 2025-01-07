
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
            # calculate LD = r^2
            snpspan = snp+1:nmarker
            ldls = cal_ldls(dosegeno, snp; snpls=snpspan)
            b =  ldls .>= minldsave               
            snpls = snpspan[b]           
            isempty(snpls) && continue
            ldls = ldls[b]

            # calculate independence LOD score
            lodls = cal_lodls(dosegeno, snp; snpls)
            b = lodls .>= minlodsave
            n = sum(b)
            n == 0 && continue

            # save results
            res = Matrix(undef,n,4)
            res[:,1] .= snp
            res[:,2] .= view(snpls,b)
            res[:,3] .= round.(view(ldls,b),digits=4)
            res[:,4] .= round.(view(lodls,b),digits=4)
            writedlm(io,res,",")
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

# calculate ld between snp and each of snpls
function cal_ldls(dosegeno::AbstractMatrix, snp::Integer;
    snpls::AbstractVector = snp+1:size(dosegeno,2))
    # 0.0 denotes missing
    nonmiss,geno_snp = findnz(dosegeno[:,snp])
    isempty(nonmiss) && return zeros(length(snpls))
    default0 = zero(eltype(geno_snp))
    cor2ls = [begin 
        nonmiss2,geno_snp2 = findnz(dosegeno[nonmiss,j])
        if !isempty(nonmiss2)            
            c = cor(geno_snp[nonmiss2], geno_snp2)
            isnan(c) ? default0 : c^2
        else
            default0
        end
    end for j in snpls]
    cor2ls
end

# calculate lod between snp and each of snpls
function cal_lodls(dosegeno::AbstractMatrix, snp::Integer;
    snpls::AbstractVector = snp+1:size(dosegeno,2), 
    doseerror::Real=0.005)    
    gstatls = cal_crosstab_gstat(dosegeno, snp; snpls,doseerror)
    lodls = [cal_linklod(i...) for i in gstatls]    
    lodls
end

# calculate gstat between snp and each of snpls
function cal_crosstab_gstat(dosegeno::AbstractMatrix, snp::Integer;
    snpls::AbstractVector = snp+1:size(dosegeno,2), 
    doseerror::Real=0.005)        
    nonmiss, geno_snp = findnz(dosegeno[:,snp])
    #  grouping columns    
    rowname = sort(unique(geno_snp))
    length(rowname)<=1 && return [zeros(2) for _ in snpls]
    rowpos = [geno_snp .== i for i in rowname]
    # remove genotypes at i with frequency <= nerror
    nerror = min(max(1,length(geno_snp)*doseerror),5)
    b = length.(rowpos) .> nerror
    rowname = rowname[b]
    rowpos = rowpos[b]
    length(rowname)<=1 && return [zeros(2) for _ in snpls]    
    res = [begin
        gcol = dosegeno[nonmiss,c]        
        colname = sort(unique(nonzeros(gcol)))
        if length(colname) <= 1
            zeros(2)
        else
            cc = [[sum(gcol[pos] .== j) for j = colname] for pos = rowpos]
            b = sum(cc) .> nerror
            cc2 = hcat(cc...)'
            cal_gstat(cc2[:,b])
        end        
    end for c in snpls]
    res
end


function cal_gstat(counts::AbstractMatrix)
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

function cal_linklod(gstat::Real,df::Real)
    df ≈ 0.0 && return 0.0
    lcdf = logcdf(Chisq(df),gstat)
    lod = invlogcdf(Chisq(1),lcdf)
    # when g is too large, lod could be NaN or Inf
    (isinf(lod) || isnan(lod)) && (lod = gstat)
    lod *= log10(ℯ)/2
    lod
end
   

# function calcrosstab(data::AbstractMatrix;doseerror::Real=0.005)
#     #  grouping rows    
#     rowname = sort(unique(data[:,1]))
#     isempty(rowname) && return [zeros(Int,0,0) for i=2:size(data,2)]
#     rowpos = [findall(data[:,1] .== i) for i= rowname]
#     # remove genotypes at i with frequency <= nerror
#     nerror = min(max(1,size(data,1)*doseerror),5)
#     pos = length.(rowpos) .> nerror
#     rowname = rowname[pos]
#     rowpos = rowpos[pos]
#     isempty(rowname) && return [zeros(Int,0,0) for i=2:size(data,2)]
#     res = [begin
#         row = view(data,:,i)
#         colname = sort(collect(skipmissing(unique(row))))
#         if isempty(colname)
#             zeros(Int,0,0)
#         else
#             cc = [[begin
#                     ls = skipmissing(row[pos])
#                     isempty(ls) ? 0 : sum(ls .== j)
#                 end for j=colname] for pos = rowpos]
#             b = sum(cc) .> nerror
#             cc2 = hcat(cc...)'
#             cc2[:,b]
#         end
#     end for i=2:size(data,2)]
#     res
# end
