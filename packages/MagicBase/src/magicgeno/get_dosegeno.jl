
function get_dosegeno(magicgeno::MagicGeno;    
    seqerror::Real=0.001,
    callthreshold::Real=0.95,
    iscalling::Bool=false,
    isdepmodel::Bool=false)
    offgeno = reduce(vcat, magicgeno.offspringgeno)
    formatls=  reduce(vcat, [i[!,:offspringformat] for i in magicgeno.markermap])
    if iscalling
        noff = size(offgeno,2)
        isoffspringinbred = isdepmodel ? trues(noff) : falses(noff)
        get_calledgeno!(offgeno,formatls; seqerror,callthreshold,isoffspringinbred)
        formatls = ["GT_unphased" for _ in 1:length(formatls)]
    end
    get_dosegeno(offgeno, formatls; isdepmodel)
end

function get_calledgeno!(calledgeno::AbstractMatrix, formatls::AbstractVector;
    seqerror::Real=0.001,
    callthreshold::Real=0.95,
    isoffspringinbred::BitVector)
    formatset = unique(formatls)
    if "GT_unphased" in formatset
        if any(isoffspringinbred)
            b = formatls .== "GT_unphased"
            subgeno =view(calledgeno, b, isoffspringinbred)
            subgeno[subgeno .== "12"] == "NN"
            subgeno[subgeno .== "21"] == "NN"
        end
    end
    if formatset != ["GT_unphased"]
        if any(isoffspringinbred)
            b = [i in ["AD","GP"] for i in formatls]
            subgeno =view(calledgeno, b, isoffspringinbred)
            subgeno .= MagicBase.genocallhaplo(subgeno, formatls[b]; seqerror,callthreshold) .^ 2
        end
        if !all(isoffspringinbred)
            b = [i in ["AD","GP"] for i in formatls]
            subgeno =view(calledgeno, b, .!isoffspringinbred)
            subgeno .= MagicBase.genocalldiplo(subgeno, formatls[b];
                seqerror = 0.001,callthreshold = 0.95)
        end
    end
    calledgeno
end


function get_dosegeno(genomtx::AbstractMatrix, formatvec::AbstractVector;
    isdepmodel::Bool=false)
    formatset = unique(formatvec)    
    # for some reasons, Int sparse matrix is slow
    float = Float64
    res = zeros(Float64,size(genomtx)...)
    for format = formatset
        if format in  ["GT_haplo", "GT_unphased","GT_phased"]
            jj = formatvec .== format
            subgeno = view(genomtx,jj,:)
            subres = view(res,jj,:)
            if !isempty(jj)
                gset = unique(subgeno)
                if isdepmodel
                    # rule = Dict("NN"=>missing,"11"=>0.0,"12"=>missing,"21"=>missing, "22"=>1.0,
                    #     "1N"=>0.0,"N1"=>0.0,"2N"=>1.0,"N2"=>1.0)                    
                    rule = Dict("11"=>1.0,"1N"=>1.0,"N1"=>1.0,"22"=>2.0,"2N"=>2.0,"N2"=>2.0)             
                    setdiff!(gset,["NN","12","21"])
                else
                    # rule = Dict("NN"=>missing,"11"=>0.0,"12"=>1.0,"21"=>1.0, "22"=>2.0,
                    #     "1N"=>0.5,"N1"=>0.5,"2N"=>1.5,"N2"=>1.5)
                    # 1+ count of allele2; 0 denotes missing
                    rule = Dict("11"=>1.0,"1N"=>1.5,"N1"=>1.5,"12"=>2.0,"21"=>2.0, "22"=>3.0,"2N"=>2.5,"N2"=>2.5)                
                    setdiff!(gset,["NN"])
                end                
                for g in gset
                    subres[subgeno .== g] .= rule[g]
                end
            end
        elseif format == "GP"
            jj = formatvec .== format
            subgeno = view(genomtx,jj,:)
            subres = view(res,jj,:)
            l = length.(subgeno)
            b = l .== 4
            subres[b] .= [max(i[1],i[2]+i[3],i[4]) <= 0.7 ? float(0.0) : float(1.0 + i[2]+i[3] + 2*i[4]) for i in subgeno[b]]
            b = l .== 3
            subres[b] .= [max(i...) <= 0.7  ? float(0.0) : float(1.0 + i[2] + 2*i[3]) for i in subgeno[b]]
            # subres[b] .= [i == [0.25,0.5,0.25]  ? float(0.0) : float(1.0 + i[2] + 2*i[3]) for i in subgeno[b]]
        else
            error(string("unknown geno format: ",format))
        end
    end
    res    
end

