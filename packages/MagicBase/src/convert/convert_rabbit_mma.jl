

"""
    rabbitped_mma2jl(mmapedfile; kwargs...)

convert input pedfile of Mathematica-version RABBIT into pedfile for Julia-version RABBIT.

# keyword arguments

`ishomozygous::Bool = false`: if true, offspring is homozygous. 

`isfglexch::Bool = false`: if true, offspring is produced with rand parent ordering.

`workdir::AbstractString=pwd()`: working directory for input and output files.

`outfile::AbstractString = "outstem_ped.csv"`: output filename.

"""
function rabbitped_mma2jl(mmapedfile::AbstractString;    
    ishomozygous::Bool = false, 
    isfglexch::Bool = false, 
    outfile::AbstractString = "outstem_ped.csv",
    workdir::AbstractString=pwd())
    pedfile = getabsfile(workdir,mmapedfile)
    isfile(pedfile) || @error string(pedfile, " does not exist")
    pedinfo = MagicBase.readmultitable(pedfile;delim=',',commentstring="##",isordereddict=true);
    infols = collect(values(pedinfo))
    length(infols ) == 2 || @error string("wrong number ",length(infols), ",  of talbes in pedinfo")
    size(infols[1],2) >=5 || @error string("wrong number ",size(infols[1],2), ", of columns in designinfo of  pedinfo")
    designdf = infols[1][!,[2,4,5]]
    if names(designdf) != ["MemberID","MotherID","FatherID"] 
        @warn string("designinfo colnames  ",names(designdf), " are not the expected ", ["MemberID","MotherID","FatherID"])
    end
    rename!(designdf,["member","mother","father"])
    size(infols[2],2) >=3 || @error string("wrong number ",size(infols[2],2), ", of columns in offspring of  pedinfo")
    offspringinfo = infols[2][!,1:3]
    if names(offspringinfo) != ["ProgenyID","MemberID","Funnelcode"]
        @warn string("offspringinfo colnames ",names(offspringinfo), " are not the expected colnames ", ["ProgenyID","MemberID","Funnelcode"])    
    end
    fcodes = unique(offspringinfo[:,3])
    if length(fcodes) > 1
        @warn string("account for different funnel codes by isfglexch=true")
        if !isfglexch
            @info string("reset isfglexch=true")
            isfglexch = true
        end
    end
    rename!(offspringinfo,["individual","member","ishomozygous"])
    offspringinfo[!,3] .= ishomozygous     
    insertcols!(offspringinfo,4,:isfglexch=>isfglexch)    
    designinfo = df2designinfo(designdf)
    founderinfo = get_founderinfo(designinfo)
    magicped = MagicPed(designinfo, founderinfo, offspringinfo)
    savemagicped(outfile, magicped;workdir)    
end

"""
    rabbitgeno_mma2jl(mmagenofile; kwargs...)

convert input genofile of Mathematica-version RABBIT into genofile for Julia-version RABBIT. 

# keyword arguments

`workdir::AbstractString=pwd()`: working directory for input and output files.

`outfile::AbstractString = "outstem_geno.vcf.gz"`: output filename.

"""
function rabbitgeno_mma2jl(mmagenofile::AbstractString;
    outfile::AbstractString = "outstem_geno.vcf.gz",
    workdir::AbstractString=pwd())
    ext = last(splitext(mmagenofile))
    ext == ".csv" || @error string("file ext = ",ext, ", not .csv")
    inputfile = getabsfile(workdir, mmagenofile)
    genodata = readdlm(inputfile,',');
    nf = genodata[1,2]
    @info string("#founders=",nf)
    genodf =DataFrame(permutedims(string.(genodata[5:end,2:end])), string.(genodata[5:end,1]))
    mapdf = DataFrame(permutedims(genodata[2:4,2:end]),["marker","linkagegroup","poscm"])
    mapdf[!,1] .= string.(mapdf[!,1])
    mapdf[!,2] .= [i == "NA" ? missing : string(i) for i in mapdf[!,2]]
    mapdf[!,3] .= [i == "NA" ? missing : Float32(i) for i in mapdf[!,3]]
    gset = unique(Matrix(genodf[!,1:nf]))
    if issubset(gset,["1","2","N"])
        fformat = "GT_haplo"    
    elseif issubset(gset,["NN", "N1", "1N", "N2", "2N", "11", "12", "21", "22"])
        fformat = "GT_unphased"
    elseif all([sum(i .== "|")==1 for i in split.(gset,"")])
        fformat = "AD"
    else
        @error string("unknown founder genotypes: ",gset)
    end
    gset = unique(Matrix(genodf[!,nf+1:end]))
    if issubset(gset,["NN", "N1", "1N", "N2", "2N", "11", "12", "21", "22"])
        offformat = "GT_unphased"
    elseif all([sum(i .== "|")==1 for i in split.(gset,"")])
        offformat = "AD"
    elseif all([in(sum(i .== "|"),[3,4]) for i in split.(gset,"")])
        offformat = "GP"
    else
        @error string("unknown offspring genotypes: ",gset)
    end    
    middf = DataFrame(physchrom = [missing for _ in 1:size(mapdf,1)], 
        physposbp=missing, info=missing,
        founderformat=fformat,offspringformat=offformat,
        foundererror=missing,offspringerror=missing,seqerror=missing,
        allelebalancemean=missing,allelebalancedisperse=missing,alleledropout=missing)
    resdf = hcat(mapdf,middf,genodf)
    isfounderinbred = fformat == "GT_haplo"
    indls = names(genodf)
    noff = length(indls) - nf
    magicped= formmagicped(JuncDist(nfounder=nf),noff)
    MagicBase.setfounderid!(magicped,indls[1:nf])
    MagicBase.setoffspringid!(magicped,indls[nf+1:end])
    magicgeno = MagicBase.formmagicgeno!(resdf, magicped;isfounderinbred, isphysmap=false)
    outfile2 = getabsfile(workdir, outfile)
    savegenodata(outfile2,magicgeno; delim='\t',commentstring = "##",workdir)
end
