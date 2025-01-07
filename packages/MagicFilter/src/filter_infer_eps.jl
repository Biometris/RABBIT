
function infer_eps(magicgeno::MagicGeno;
    model::AbstractString,
    snpthin::Integer=1,
    mapavailable::Bool=false,
    accuracy_eps::Real = 0.001,
    maxiteration::Integer = 20,
    isfounderinbred::Bool=true,    
    historyfile::AbstractString,
    logio::Union{Nothing,IO},
    verbose::Bool)
    if snpthin > 1
        msg = string("take every ",snpthin, "-th marker")
        MagicBase.printconsole(logio,verbose, msg)
        snpsubset = 1:snpthin:max(size.(magicgeno.markermap,1)...)
        subgeno = submagicgeno(magicgeno;snpsubset)
    else
        subgeno = magicgeno
    end
    tused = @elapsed magicprior=MagicReconstruct.calmagicprior(subgeno,model;isfounderinbred)
    msg = string("tused=",round(tused,digits=2), "s in calculating prior")
    MagicBase.printconsole(logio,verbose, msg)
    tused = @elapsed pretuplels = calpretuplels(subgeno, magicprior;
        model,isfounderinbred);
    msg = string("tused=",round(tused,digits=2), "s in precompution")
    MagicBase.printconsole(logio,verbose, msg)
    nfounder = size(subgeno.magicped.founderinfo,1)
    ind2subpop =  subgeno.magicped.offspringinfo[!,:member]
    noff = length(ind2subpop)
    epsf = epso = 0.005
    epsf_vec = epsf*ones(nfounder)
    epso_vec = epso*ones(noff)
    oldepsf = zeros(nfounder)
    oldepso = zeros(noff)
    loglls = zeros(noff)
    open(historyfile,"w") do his_io
        write(his_io, "founderlist: ", join(subgeno.magicped.founderinfo[!,:individual],","),"\n")
        write(his_io, "offspringlist: ", join(subgeno.magicped.offspringinfo[!,:individual],","),"\n")
        write(his_io, "iteration, epsf_vec, epso_vec\n")
        write(his_io, "0,", join(vcat(epsf_vec,epso_vec),","),"\n")
        brentmaxiter = 10
        brentstepsize = 4.6 # exp(2.3) ~ 10, exp(1.6)~ 5
        for it in 1:maxiteration
            startt = time()
            oldepsf .= epsf_vec
            oldepso .= epso_vec
            @showprogress 1 "calculating epsf..." for findex in 1:nfounder
                # println("findex=",findex)
                update_epsf_findex!(epsf_vec,epso_vec,findex,
                    pretuplels,mapavailable,brentmaxiter,brentstepsize)
            end
            progress = Progress(noff;desc="calculating epso...")       
            ThreadsX.foreach(eachindex(epso_vec, loglls, ind2subpop)) do off                 
                epso_vec[off], loglls[off] = infer_epso_offspring(epsf_vec,
                    epso_vec[off],off,ind2subpop[off],
                    pretuplels,mapavailable, brentmaxiter,brentstepsize)
                # println("off=",off, "; epso_vec[off], loglls[off] = ", [epso_vec[off], loglls[off]])
                next!(progress)
            end
            logl = sum(loglls)
            diff_epsf = max((epsf_vec .- oldepsf)...)
            diff_epso = max((epso_vec .- oldepso)...)
            diff_eps = round.([diff_epsf,diff_epso]; digits=4)
            maxeps = round.([max(epsf_vec...),max(epso_vec...)],digits=2)
            msg = string("it=",it,
                ", diff_eps=", join(diff_eps,"|"),
                ", max_eps=", join(maxeps,"|"),
                ", logl=",round(logl,digits=1),
                ", tused=",round(time()-startt,digits=1),"s")
            MagicBase.printconsole(logio,verbose, msg)
            write(his_io, string(it), ",", join(vcat(epsf_vec,epso_vec),","),"\n")
            max(diff_epsf,diff_epso) < accuracy_eps && break
        end
    end
    epsf_vec,epso_vec
end

function infer_epso_offspring(epsf_vec::AbstractVector,epso::Real,
    offspring::Integer,subpopid::AbstractString,
    pretuplels::AbstractVector, mapavailable::Bool,
    maxiter::Integer,stepsize::Real)
    # epso refers to the offspring
    accuracygoal, precisiongoal = 2, 2
    xstart = log(max(10^(-5),epso))
    lowbound = max(log(10^(-5)),xstart-stepsize)
    upbound = min(0,xstart+stepsize)
    function loglfun(x::Real)
        loglike_ind(epsf_vec,exp(x),offspring,subpopid,pretuplels,mapavailable)
    end
    res= MagicBase.brentMax(loglfun,lowbound,upbound;
        xstart, precisiongoal,accuracygoal,maxiter)
    newepso, logl = exp(res[1]), res[2]
    newepso, logl
end

function update_epsf_findex!(epsf_vec::AbstractVector,epso_vec::AbstractVector,
    findex::Integer,pretuplels::AbstractVector,mapavailable::Bool,
    maxiter::Integer,stepsize::Real)
    accuracygoal, precisiongoal = 2, 2
    xstart = log(max(1e-5,epsf_vec[findex]))
    lowbound = max(log(1e-5),xstart-stepsize)
    upbound = min(0,xstart+stepsize)
    popidls = popfromfindex(findex,first(pretuplels).popmakeup)
    function loglfun(x::Real)
        epsf_vec[findex] = exp(x)
        loglike_popls(epsf_vec,epso_vec,popidls,pretuplels,mapavailable)
    end
    res= MagicBase.brentMax(loglfun,lowbound, upbound;
        xstart, precisiongoal,accuracygoal,maxiter)
    epsf_vec[findex] = exp(res[1])
    res[2]
end

function popfromfindex(findex::Integer,popmakeup::AbstractDict)
    popfromfindex([findex],popmakeup)
end

# list of all subpopulations whose parnets are contained in the vector findex
function popfromfindex(findex::AbstractVector,popmakeup::AbstractDict)
    res = []
    for (key, val) =popmakeup
        a = intersect(findex,val["founder"])
        isempty(a) || push!(res,key)
    end
    res
end
