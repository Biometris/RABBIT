
struct FglGeno
    fgl::Vector
    genetic_value::Float64
end

function setfounderfgl(ped::Pedigrees.Pedigree,chrlen::AbstractVector,
    isfounderinbred::Bool,nodesize::Integer)
    nmem = length(ped.member)
    pedfgl = Vector{Vector{FglGeno}}(undef,nmem)
    nf = ped.nfounder
    @assert eltype(ped.member)<: Integer
    for i=1:nf
        mem = ped.member[i]
        fgl = isfounderinbred ? [mem,mem] : [2mem-1,2mem]
        fgl2 = [[[[0.0,fgl[o]], [len,NaN]] for o=1:2] for len=chrlen]
        genetic_value = 0.0
        pedfgl[i] = [FglGeno(fgl2,genetic_value) for _ in 1:nodesize]
    end
    pedfgl
end

function randbreaks(chrlen::Real,isobligate::Bool, interference::Integer)
    nbreak = (1+interference)*rand(Poisson(chrlen/100))
    isobligate && (nbreak = max(1, nbreak))
    breaks = sort(rand(Uniform(0,chrlen),nbreak))
    if length(breaks) >= 3
        breaks = breaks[rand(1:1 + interference):1+interference:end]
    end
    pushfirst!(breaks,0)
    push!(breaks,chrlen)
    breaks
end

function breakhomolog(homolog::AbstractVector,breaks::AbstractVector)
    nbreak = length(breaks)
    # res_cols: breakpoint, fgl, segment
    res = Matrix(undef,nbreak,3)
    res[:,1]=breaks
    res[1,2] = homolog[1][2]
    res[end,2] = homolog[end][2]
    res[end,3] = []
    h=2
    for i = 1:nbreak-1
        res[i, 3] = []
        while h < length(homolog)
            if homolog[h][1] < res[i+1, 1]
                push!(res[i, 3], homolog[h])
                h += 1
            else
                break
            end
        end
    end
    for i= 1:nbreak-2
        res[i+1,2] = isempty(res[i,3]) ? res[i,2] : res[i,3][end][2]
    end
    res
end

# maternal = [[0.0, 2.0], [200.0, -1.0]]
# paternal = [[0.0, 3.0],[30.9,2.0], [128.4, 4], [200.0, -1.0]]
# homologpair = [maternal,paternal]
# gamete = crossover(homologpair, isobligate,interference)
# hcat(gamete...)
function crossover(homologpair::AbstractVector,isobligate::Bool,interference::Integer)
    maternal, paternal = homologpair
    chrlen = maternal[end][1]
    # @assert chrlen ≈ paternal[end][1]
    breaks = randbreaks(chrlen, isobligate,interference)
    breakm = breakhomolog(maternal,breaks)
    breakp = breakhomolog(paternal,breaks)
    oo = trues(length(breaks))
    i0 = rand()<0.5 ? 1 : 2
    oo[i0:2:end] .= false
    res=Vector{Vector{Float64}}()
    for i in eachindex(oo)
        if oo[i]
            push!(res,[breakm[i,1],breakm[i,2]])
            push!(res,breakm[i,3]...)
        else
            push!(res,[breakp[i,1],breakp[i,2]])
            push!(res,breakp[i,3]...)
        end
    end
    # writedlm(stdout,hcat(res...))
    # show(stdout, "text/plain", hcat(res...))
    for i=2:length(res)-1
        res[i][2]==res[i-1][2] && (res[i] = res[i-1])
    end
    unique(res)
end

function simpedfgl(ped::Pedigrees.Pedigree,chrlen::AbstractVector;
    isfounderinbred::Bool,
    isobligate::Bool,
    interference::Integer,
    select_trait::Union{Nothing,TraitQtl},
    select_prop::Real,
    nodesize_aft::Integer)
    # gender must be  0=notapplicable, 1=female or 2=male
    # If isoogamy=True,the last pair of chromosomes are sex chromosomes.
    # XY for male and XX for female.
    # linkage group = [maternally derived chromosome, paternally derived chromosome]
    isoogamy = (2 in unique(ped.gender) && length(chrlen)>1)
    pedfgl = setfounderfgl(ped,chrlen,isfounderinbred,nodesize_aft)
    nf = ped.nfounder
    memberindices = zeros(Int,max(ped.member...))
    memberindices[ped.member] = 1:length(ped.member)
    nodesize_bef = round(Int, nodesize_aft/select_prop)
    for ind = nf+1:length(pedfgl)
        motherindex = memberindices[ped.mother[ind]]
        fatherindex = memberindices[ped.father[ind]]
        fgl_dam = pedfgl[motherindex]
        fgl_sire = pedfgl[fatherindex]
        @assert length(fgl_dam) == length(fgl_sire) == nodesize_aft
        if motherindex == fatherindex
            indexls = [(i,i) for i in 1:nodesize_aft]
        else
            if nodesize_bef > nodesize_aft^2
                indexls = [(rand(1:nodesize_aft),rand(1:nodesize_aft)) for _ in 1:nodesize_bef]
            else
                indexls = sample(collect(Iterators.product(1:nodesize_aft,1:nodesize_aft)),
                    nodesize_bef;replace=false)
            end
        end
        fglnode = Vector{FglGeno}(undef,length(indexls))
        for i in eachindex(indexls)
            idam,isire = indexls[i]
            egg  = crossover.(fgl_dam[idam].fgl, isobligate,interference)
            sperm  = crossover.(fgl_sire[isire].fgl, isobligate,interference)
            zygote = map((x,y)->[x,y],egg,sperm)
            if isoogamy
                # if isoogamy, gender 1=female or 2=male; 0=notapplicable is not allowed
                # sex chromsome is in the last
                gender = ped.gender[ind]
                sexsperm = pedfgl[father[1]][end][gender]
                zygote[end][2] = sexsperm
            end
            genetic_value = select_prop==1.0 ? 0.0 : trait_genetic_val(zygote,select_trait)
            fglnode[i] = FglGeno(zygote,genetic_value)
        end
        oo = sortperm([i.genetic_value for i in fglnode],rev=true)
        pedfgl[ind] = fglnode[oo[1:nodesize_aft]]
    end
    pedfgl
end

function getpeddict(magicped::MagicPed)
    designped = magicped.designinfo
    indexped = Pedigrees.toindexped(designped)
    memberdict = Dict(designped.member .=> indexped.member)
    memberls  = magicped.offspringinfo[!,:member]
    subpopls = unique(memberls)
    @assert issubset(subpopls,designped.member)
    peddict = Dict([begin
        offls = findall(memberls .== subpop)
        subpop2 = get(memberdict,subpop,nothing)
        subped = Pedigrees.getsubped(indexped,subpop2)
        subpop => (offls,subped)
    end for subpop = subpopls])
    peddict
end

function simcontfgl(magicped::MagicPed;
    chrid::AbstractVector=["chr1"],
    chrlen::AbstractVector=[100], # centiMorgan
    isfounderinbred::Bool=true,
    isobligate::Bool=false,
    interference::Integer=0,
    select_nqtl::Integer=50,
    select_dom::Real=0, 
    select_prop::Real = 1.0)
    @assert length(chrid) == length(chrlen)
    designped = magicped.designinfo
    indexped = Pedigrees.toindexped(designped)
    nf = indexped.nfounder
    if select_prop == 1.0 
        select_trait =  nothing 
    else
        select_trait = TraitQtl(nf, chrlen, select_nqtl, select_dom)
    end    
    nondiv = allunique(magicped.offspringinfo[!,:member])
    if nondiv
        pedfgl =simpedfgl(indexped,chrlen;
            isfounderinbred,isobligate,interference,
            select_trait, select_prop,nodesize_aft=1)
        founderfgl = [first(i).fgl for i in pedfgl[1:nf]]
        memberdict = Dict(designped.member .=> indexped.member)
        membls =  [get(memberdict,i,nothing) for i=magicped.offspringinfo[!,:member]]
        offspringfgl = [rand(i).fgl for i in pedfgl[membls]]
    else
        tempfgl = setfounderfgl(indexped,chrlen,isfounderinbred,1)
        founderfgl = [first(i).fgl for i in tempfgl[1:nf]]
        noff  = size(magicped.offspringinfo,1)
        offspringfgl = Vector(undef,noff)
        peddict = getpeddict(magicped)
        for subpop = keys(peddict)
            offls, subped = peddict[subpop]
            subnoff = length(offls)
            lastnode = last(simpedfgl(subped,chrlen;
                isfounderinbred,isobligate,interference,
                select_trait, select_prop,nodesize_aft=subnoff))
            offspringfgl[offls] = [i.fgl for i in lastnode]
        end
    end
    markermap = MagicBase.contfgl_markermap(chrlen,chrid)
    nchr = length(chrlen)    
    # foundergeno = MagicBase.contfgl_foundergeno(chrlen,isfounderinbred,nf)
    foundergeno  = [hcat([i[chr] for i=founderfgl]...) for chr=1:nchr]
    offspringgeno  = [hcat([i[chr] for i=offspringfgl]...) for chr=1:nchr]
    misc = Dict{String, DataFrame}()
    contfgl  = MagicGeno(magicped,markermap,foundergeno,offspringgeno,misc)
    simfgl_DH!(contfgl)
    contfgl
end

function simfgl_DH!(contfgl::MagicGeno)
    homooffls = findall(contfgl.magicped.offspringinfo[!,:ishomozygous])
    for offgeno in contfgl.offspringgeno
        for off in homooffls
            if rand()<= 0.5
                offgeno[1,off] = copy(offgeno[2,off])
            else
                offgeno[2,off] = copy(offgeno[1,off])
            end
        end
    end
end

#############################################################

function gridhomologfgl(homologfgl::Vector{Vector{T}} where T <: AbstractFloat,
    chrmarkerpos::AbstractVector)
    # @assert chrmarkerpos[end] < homologfgl[end][1]
    res = zeros(Int,length(chrmarkerpos))
    seg = 1
    i=1
    while i<=length(chrmarkerpos)
        if chrmarkerpos[i]<homologfgl[seg+1][1]
            res[i] = homologfgl[seg][2]
            i += 1
        else
            seg += 1
        end
    end
    # pos = findall(abs.(diff(res)) .> 0)
    # [chrmarkerpos[pos] chrmarkerpos[pos .+ 1]]
    # scatter(map((x,y)->(x,y), chrmarkerpos, res))
    res
end

function getgridfgl(offspringfgl::Vector{Matrix},markerpos::AbstractVector)
    nchr= length(markerpos)
    @assert length(offspringfgl) == nchr
    gridfgl  = [begin
        chrfgl = [gridhomologfgl(i,markerpos[chr]) for i = offspringfgl[chr]]
        # @assert size(chrfgl,1)==2
        chrfgl2 = [let (hm, hp) = col
            [[hm[i], hp[i]] for i = 1:length(hm)]
        end for col in eachcol(chrfgl)]
        hcat(chrfgl2...)
    end for chr=1:nchr]
    gridfgl
end

function checkfounder!(founderhaplo::MagicGeno, magicped::MagicPed)
    # modeify founderhaplo for consistent parent IDs and ordering
    parents = magicped.founderinfo[!,:individual]
    parents2 = founderhaplo.magicped.founderinfo[!,:individual]
    diff = setdiff(parents,parents2)
    if !isempty(diff)
        msg = string("parents in magicped but not in founderhaplo:",diff)
        error(msg)
    end
    fdict = Dict(parents2 .=> 1:length(parents2))
    indices = [get(fdict,i,nothing) for i=parents]
    @assert parents2[indices] == parents
    if indices == 1:length(parents)
        founderhaplo.magicped.founderinfo = founderhaplo.magicped.founderinfo[indices,:]
        founderhaplo.foundergeno = [geno[:,indices] for geno =founderhaplo.foundergeno]
    end
    founderhaplo
end


function simfgl(founderhaplo::MagicGeno,magicped::MagicPed;
    isfounderinbred::Bool=true,
    isobligate::Bool=false,
    interference::Integer=0,
    select_nqtl::Integer=50,    
    select_prop::Real = 0.1)
    checkfounder!(founderhaplo,magicped)
    markermap = deepcopy(founderhaplo.markermap)
    # chrlen+1: ensure chrlen is larger than the last marker position
    chrid = [i[1,:linkagegroup] for i=markermap]
    chrlen = [ceil(i[end,:poscm]-i[1,:poscm])+1 for i=markermap]
    markerpos = [i[!,:poscm] .- i[1,:poscm] for i=markermap]
    contfgl = simcontfgl(magicped; chrid,chrlen,isfounderinbred,isobligate,interference,
        select_nqtl,select_prop)    
    offgeno = getgridfgl(contfgl.offspringgeno, markerpos)
    # set fgl definition in mics
    fgl = Int.([i[1][2] for i=contfgl.foundergeno[1]])
    founderid = contfgl.magicped.founderinfo[!,:individual]
    df = DataFrame(:founder =>founderid, :maternal=>fgl[1,:],:paternal =>fgl[2,:])
    misc = Dict("fgl"=>df)
    # set format
    for i=markermap
        i[!,:offspringformat] .= "discretefgl"
    end
    magicfgl = MagicGeno(magicped,markermap,founderhaplo.foundergeno,offgeno,misc)
    contfgl, magicfgl
end
