# diploprior returns P(Z|D, O, epsf) at a given site,
#  where Z=true genotype, D=derived genotype,
#  O= ancestral origins, epsf = founder allelic error probability
# return prior with size: 4 true genotypes x 9 derived genotypes
function diploprior(epsf::Real)
    # 9 derived genotypes =Z: "NN", "N1", "1N", "N2", "2N", "11", "12", "21", "22"
    # 4 true genotypes=Z: "11", "12", "21", "22"
    ptrueibd = zeros(9,4)
    ptruenonibd = zeros(9,4)
    ptrueibd[1,:] .= [1, 0, 0, 1]
    ptrueibd[6,:] .= [1-epsf,0,0,epsf]
    ptrueibd[9,:] .= [epsf,0,0,1-epsf]
    for i=[1,6,9]
        ptrueibd[i,:] ./= sum(ptrueibd[i,:])
    end
    ptruenonibd[1,:] .= [1, 1, 1, 1]
    ptruenonibd[2,:] .= [1-epsf, epsf, 1-epsf, epsf]
    ptruenonibd[3,:] .= [1-epsf, 1-epsf, epsf, epsf]
    ptruenonibd[4,:] .= [epsf, 1-epsf, epsf, 1-epsf]
    ptruenonibd[5,:] .= [epsf, epsf, 1-epsf, 1-epsf]
    ptruenonibd[6,:] .= [(1-epsf)^2, (1-epsf)epsf, (1-epsf)epsf, epsf^2]
    ptruenonibd[7,:] .= [(1-epsf)epsf, (1-epsf)^2, epsf^2, (1-epsf)epsf]
    ptruenonibd[8,:] .= [(1-epsf)epsf, epsf^2, (1-epsf)^2, (1-epsf)epsf]
    ptruenonibd[9,:] .= [epsf^2, (1-epsf)epsf, (1-epsf)epsf, (1-epsf)^2]
    for i=1:9
        ptruenonibd[i,:] ./= sum(ptruenonibd[i,:])
    end
    (nonibd=permutedims(ptruenonibd),ibd=permutedims(ptrueibd))
end

function diploprior(origin_pair::AbstractVector, epsf_pair::AbstractVector)
    # 9 derived genotypes =Z: "NN", "N1", "1N", "N2", "2N", "11", "12", "21", "22"
    # 4 true genotypes=Z: "11", "12", "21", "22"
    ptrue = zeros(9,4)
    isibd = isequal(origin_pair...)
    if isibd
        isequal(epsf_pair...) || @error string("epsf_pair=",epsf_pair," not equal for ibd state")
        epsf = first(epsf_pair)
        ptrue[1,:] .= [1, 0, 0, 1]
        ptrue[6,:] .= [1-epsf,0,0,epsf]
        ptrue[9,:] .= [epsf,0,0,1-epsf]
        for i in [1,6,9]
            ptrue[i,:] ./= sum(ptrue[i,:])
        end
    else
        e1,e2 = epsf_pair
        ptrue[1,:] .= [1, 1, 1, 1]
        ptrue[2,:] .= [1-e2, e2, 1-e2, e2]
        ptrue[3,:] .= [1-e1, 1-e1, e1, e1]
        ptrue[4,:] .= [e2, 1-e2, e2, 1-e2]
        ptrue[5,:] .= [e1, e1, 1-e1, 1-e1]
        ptrue[6,:] .= [(1-e1)*(1-e2), (1-e1)*e2, e1*(1-e2), e1*e2]
        ptrue[7,:] .= [(1-e1)*e2, (1-e1)*(1-e2), e1*e2, e1*(1-e2)]
        ptrue[8,:] .= [e1*(1-e2), e1*e2, (1-e1)*(1-e2), (1-e1)*e2]
        ptrue[9,:] .= [e1*e2, e1*(1-e2), (1-e1)*e2, (1-e1)*(1-e2)]
        for i in 1:9
            ptrue[i,:] ./= sum(ptrue[i,:])
        end
    end
    ptrue'
end

function diplolike(epso::Real; isoffphased::Bool=false, israndallele::Bool)    
    if israndallele
        # genotyping error model: random allelic model
        if isoffphased
            # 9 observed phased genotypes =Y: "NN", "N1", "1N", "N2", "2N", "11", "12", "21", "22"
            # 4 true genotypes=Z: "11", "12", "21", "22"
            like = zeros(9,4)
            like[1,:] .= [1, 1, 1, 1]
            like[2,:] .= [1-epso, 0.5, 0.5, epso]
            like[3,:] .= [1-epso, 0.5, 0.5, epso]
            like[4,:] .= [epso, 0.5, 0.5, 1-epso]
            like[5,:] .= [epso, 0.5, 0.5, 1-epso]
            like[6,:] .= [(1-epso)^2, (1-epso)*epso, (1-epso)*epso, epso^2]
            like[7,:] .= [(1-epso)*epso, (1-epso)^2, epso^2, (1-epso)*epso]
            like[8,:] .= [(1-epso)*epso, epso^2, (1-epso)^2, (1-epso)*epso]
            like[9,:] .= [epso^2, (1-epso)*epso, (1-epso)*epso, (1-epso)^2]
        else
            # 6 observed unphased genotypes =Y: "NN", "N1", "N2", 11", "12", "22"
            # 4 true genotypes=Z: "11", "12", "21", "22"
            like = zeros(6,4)
            like[1,:] .= [1, 1, 1, 1]
            like[2,:] .= [1-epso, 0.5, 0.5, epso]
            like[3,:] .= [epso, 0.5, 0.5, 1-epso]
            like[4,:] .= [(1-epso)^2, (1-epso)*epso, (1-epso)*epso, epso^2]
            like[5,:] .= [2*(1-epso)*epso, (1-epso)^2+epso^2, (1-epso)^2+epso^2, 2*(1-epso)*epso]
            like[6,:] .= [epso^2, (1-epso)*epso, (1-epso)*epso, (1-epso)^2]
        end
    else
        # genotyping error model: random genotypic model
        if isoffphased
            # 9 observed phased genotypes =Y: "NN", "N1", "1N", "N2", "2N", "11", "12", "21", "22"
            # 4 true genotypes=Z: "11", "12", "21", "22"
            like = zeros(9,4)
            like[1,:] .= [1, 1, 1, 1]
            like[2,:] .= [1-epso, 0.5, 0.5, epso]
            like[3,:] .= [1-epso, 0.5, 0.5, epso]
            like[4,:] .= [epso, 0.5, 0.5, 1-epso]
            like[5,:] .= [epso, 0.5, 0.5, 1-epso]
            like[6,:] .= [1-epso, epso/3, epso/3, epso/3]
            like[7,:] .= [epso/3, 1-epso, epso/3, epso/3]
            like[8,:] .= [epso/3, epso/3, 1-epso, epso/3]
            like[9,:] .= [epso/3, epso/3, epso/3, 1-epso]
        else
            # 6 observed unphased genotypes =Y: "NN", "N1", "N2", 11", "12", "22"
            # 4 true genotypes=Z: "11", "12", "21", "22"
            like = zeros(6,4)
            like[1,:] .= [1, 1, 1, 1]
            like[2,:] .= [1-epso, 0.5, 0.5, epso]
            like[3,:] .= [epso, 0.5, 0.5, 1-epso]
            like[4,:] .= [1-epso, epso/3, epso/3, epso/3]
            like[5,:] .= [epso*2/3, 1-epso*2/3, 1-epso*2/3, epso*2/3]
            like[6,:] .= [epso/3,epso/3,epso/3,1-epso]
        end   
    end
    like
end

# return like with size 6 (or 9 for phased) obsvered genotypes x 9 derived genotypes
function diplolike(epsf::Real, epso::Real; isoffphased::Bool=false, israndallele::Bool)    
    prior=diploprior(epsf)
    like=diplolike(epso; isoffphased,israndallele)
    (nonibd=like*prior.nonibd, ibd=like*prior.ibd)
end

# 3 derived haplotypes: N, 1, 2. They are NN, 11, 22 for inbred populations
# 2 true haplotypes: 1, 2. They are 11, 22 for inbred populations
# return prior with size 2 true haplotypes x 3 derived haplotypes
function haploprior(epsf::Real)
    [0.5 1.0-epsf epsf; 
     0.5 epsf 1.0-epsf]
end

# return like with size 3 obsvered haplotypes x 2 true haplotypes
haplolike(epso::Real) = [1.0 1.0; 1.0-epso epso; epso 1.0-epso]

# return like with size 3 obsvered haplotypes x 3 derived haplotypes
function haplolike(epsf::Real, epso::Real)
    prior=haploprior(epsf)
    like=haplolike(epso)
    like*prior
end

function calfderive(fhaplo::AbstractVector;
    nzstate::Union{Nothing,AbstractVector}=nothing,
    ishaploid::Bool=false)
    vec(calfderive(reshape(fhaplo,1,:);nzstate,ishaploid))
end

# fhaplo[m,p] = founder allele at marker m for inbredparent p.
# or fhaplo[phase,p] = founder allele for inbredparent p for around pahses at a given site
# fderive size: (nsnp, nstate), where nstate = size(fhaplo,2) for ishaploid=true, and size(fhaplo,2)^2 for ishaploid =false
function calfderive(fhaplo::AbstractMatrix;
    nzstate::Union{Nothing,AbstractVector}=nothing,
    ishaploid::Bool=false)
    ty=typeof(fhaplo[1,1])
    if ty <: AbstractString
        # unphased founder genotypes, eg, "12", or "1", "2"
        fhaplo2 = fhaplo
    elseif ty <: AbstractVector
        # phased founder genotypes, eg,  ["1", "2"]
        fhaplo2 = permutedims(reduce(vcat, [reduce(hcat,i) for i=eachcol(fhaplo)]))
    else
        @error string("unknown chrfhaplo eltype: ",ty)
    end
    if ishaploid
        derrule=Dict(["N"=>Int8(1),"1"=>Int8(2), "2"=>Int8(3)])
        fderive=[get(derrule,i, Int8(-1)) for i = fhaplo2]
        if -1 in unique(fderive)
            error("there exist unknow founder alleles")
        end
    else
        # unphased founder genotypes, eg, "12"                
        nsnp,nfgl=size(fhaplo2)
        state = MagicBase.prior_diploindex(nfgl)
        isnothing(nzstate) && (nzstate = 1:length(state))
        state = state[nzstate]
        # try 
        #     state = state[nzstate]
        # catch err
        #     println("nfgl=",nfgl,",state=",state,",nzstate=",nzstate)
        #     @info string(err)
        # end
        derrule=Dict(["NN", "N1", "1N", "N2", "2N", "11", "12", "21", "22"] .=> Int8.(1:9))
        # fderive=Matrix{Union{Missing,Int}}(missing,nsnp,nfgl^2)
        fderive=spzeros(Int8,nsnp,nfgl^2)        
        fderive[:,nzstate]=[get(derrule,join(fhaplo2[snp,i]),Int8(-1)) for snp=1:nsnp,i=state]
        if -1 in unique(fderive[:,nzstate])
            error("there exist unknown founder alleles")
        end
    end
    fderive
end

# obsgeno[i] = observed genotypes of offpsrig i
# all offspring have the same ishaploid
function caloffcode(obsgeno::AbstractVector;
    ishaploid::Bool=false,isoffphased::Bool=false)
    obsgeno2 = reshape(obsgeno,1,:)
    offcode = caloffcode(obsgeno2; ishaploid,isoffphased)
    offcode[1,:]
end

# obsgeno[t, i] = observed genotypes of offpsrig i at locus t
# all offspring (columns) have the same is haploid
function caloffcode(obsgeno::AbstractMatrix;
    ishaploid::Bool=false,isoffphased::Bool=false)    
    if ishaploid
        dict=Dict(["N","1","2"] .=> Int8.(1:3))        
        merge!(dict,Dict(["NY","1Y","2Y"] .=> Int8.(1:3)))
        merge!(dict,Dict(["NN","12","21","1N","N1","11","N2","2N","22"] .=> Int8.([1,1,1,2,2,2,3,3,3])))
        offcode=[get(dict,isoffphased ? join(sort(i)) : i,missing) for i=obsgeno]
        b = ismissing.(offcode)        
        if any(b)
            error("unknown observed offspringgeno: ", unique(obsgeno[b]))
        end
    else
        if isoffphased
            dict=Dict([["N","N"], ["1","N"], ["N","1"], ["2","N"], ["N","2"],
                ["1","1"], ["1","2"],["2","1"], ["2","2"]] .=> Int8.(1:9))
            merge!(dict,Dict([["N","Y"]=>Int8(1),["1","Y"]=>Int8(6),["2","Y"]=>Int8(9)]))
        else
            dict = Dict(["NN", "1N", "2N", "11", "12", "22"] .=> Int8.(1:6))
            merge!(dict,Dict(["N1"=>Int8(2),"N2"=>Int8(3),"21"=>Int8(5)]))
        end
        offcode =[get(dict,i,missing) for i=obsgeno]        
        b = ismissing.(offcode)
        if any(b)
            error("unknown observed offspringgeno: ",unique(obsgeno[b]), ", isoffphased=",isoffphased)
        end
    end
    offcode
end

function callinelikesnp!(dataprobls::AbstractVector,
    condlike::Union{NamedTuple,AbstractMatrix},
    fderive::AbstractMatrix, nzstate::AbstractVector, obsseq::AbstractVector; 
    ishaploid::Bool=false)
    condlikels = repeat([condlike],length(obsseq))
    callinelikesnp!(dataprobls,condlikels,fderive,nzstate, obsseq; ishaploid)
end

function callinelikesnp!(dataprobls::AbstractVector,
    condlike::Vector{T} where T <: Union{NamedTuple,AbstractMatrix},
    fderive::AbstractMatrix, nzstate::AbstractVector, obsseq::AbstractVector; 
    ishaploid::Bool=false)
    length(obsseq) == length(condlike) == size(fderive,1) || @error "dimensions mismatch."
    if ishaploid
        for snp in eachindex(obsseq)
            dataprobls[snp] .= condlike[snp][obsseq[snp],fderive[snp,nzstate]] 
        end
    else        
        nfgl0=sqrt(size(fderive,2))
        if nfgl0 % 1 !=0
            error("unexpected number of columns in fderive")
        end
        nfgl=Int(nfgl0)       
        ibdls = [(i-1)nfgl+i for i in 1:nfgl]
        isnzibd = [i in ibdls for i in nzstate]
        nzibdls = nzstate[isnzibd]
        nznonibdls = nzstate[.!isnzibd]  
        # check any missing of fderive at nzstate
        # missing pattern is the same for each row
        # missing beeing zero in the sparse matrix of fderive
        all(fderive[1,nzstate] .> 0) || @error "inconsistent nzstate and fderive"
        issparse = length(nzstate)/(nfgl^2) < 0.5 
        for snp in eachindex(obsseq)
            prob = issparse ? spzeros(_float_like, nfgl^2) : zeros(_float_like,nfgl^2)
            prob[nzibdls] .= condlike[snp].ibd[obsseq[snp],fderive[snp,nzibdls]]
            prob[nznonibdls] .= condlike[snp].nonibd[obsseq[snp],fderive[snp,nznonibdls]]
            dataprobls[snp] .= Vector(prob[nzstate])
        end 
    end    
end
