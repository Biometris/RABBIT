
"""
    parsebreedped(pedfile; fixed_nself=10, outfile, commentstring="##",workdir=pwd())

convert a breedped pedfile into the magicped outfile. The first 3 columns of breedpedfile must be sample, pedcode, nself.

# Keyword arguments

`fixed_nself::Integer = 20`: interprete "FIXED" in nself column as 20. 

`delim::AbstractChar=','`: text delimitor of input pedfile. 

`commentstring::AbstractString="##"`: the lines beginning with commentstring are ignored
in pedfile.

`workdir::AbstractString=pwd()`: directory for reading pedfile.

"""
function parsebreedped(pedfile::AbstractString;
    fixed_nself::Integer = 20,
    commentstring::AbstractString="##",
    delim::AbstractChar=',',    
    outfile::Union{Nothing,AbstractString} = first(splitext(basename(pedfile)))*"_magicped.csv",
    workdir::AbstractString=pwd())
    # pedfile contains 3 columns: sample, pedcode, nself
    # S0 denotes the individual resulting from the last cross, Sn=nself
    ext = last(split_allext(pedfile))
    ext in [".csv",".csv.gz"] || @error string("pedfile = ",pedfile, 
        ", file extension = ", ext, ", is not in .csv, or .csv.gz")
    pedfile2=getabsfile(workdir,pedfile)
    isfile(pedfile2) || @error string(pedfile2, " does not exist")
    peddf = CSV.read(pedfile2,DataFrame; delim,comment=commentstring)
    pedcode = peddf[!,2]
    gencol = uppercase.(strip.(string.(peddf[!,3])))
    ishomozygous = [g in ["DH"] for g in gencol]        
    ngeration = [if g in ["DH"]
            g
        elseif  in(g,["FIXED"])
            fixed_nself
        elseif occursin(r"^[sS]([0-9]{1,})",g)
            m = match(r"^[sS]([0-9]{1,})",g)
            g2 = only(m.captures)
            parse(Int,g2)
        elseif occursin(r"^[0-9]{1,}",g)
            parse(Int,g)
        else
            @error string("Could not parse number of inbreeding generations: ",g)
        end for g in gencol]
    dfls = map((x,y)->parsebreedcode(string(x,"=>", y); fixed_nself),pedcode,ngeration)
    member = [i[end,:member] for i in dfls]
    df = unique(reduce(vcat, dfls))
    designinfo = Pedigree(df)
    nf=designinfo.nfounder
    founderinfo = DataFrame(:individual=>designinfo.member[1:nf],
        :gender=>designinfo.gender[1:nf])
    ind = string.(peddf[:,1])
    dict = Dict(designinfo.member .=> 1:length(designinfo.member))
    ii = [get(dict, i, nothing) for i in member]
    gender = designinfo.gender[ii]
    offspringinfo = DataFrame([:individual=> ind,:member=>member,
        :ishomozygous=>ishomozygous,
        :isfglexch=>false,
        :gender=>gender])
    magicped = MagicPed(designinfo,founderinfo,offspringinfo)
    if !isnothing(outfile)
        savemagicped(outfile,magicped; workdir)
    end
    magicped
end


function generate_breedped(;
    subpopsize::Integer = 10, 
    outfile::AbstractString = "breedped.csv",
    workdir::AbstractString=pwd())
    sampleid = []
    pedcode = []
    nself = []
    # F2: P1, P2
    append!(pedcode,["P1/4/P1/3/P1//P1/P2" for _ in 1:subpopsize])
    append!(nself,[2 for _ in 1:subpopsize])
    append!(sampleid,[string("BC_",i) for i in 1:subpopsize])
    # DH: P2, P3
    append!(pedcode,["P2/P3" for _ in 1:subpopsize])
    append!(nself,["DH" for _ in 1:subpopsize])
    append!(sampleid,[string("DH_",i) for i in 1:subpopsize])
    # 4ril-self4: P3, P4, P5, P6
    append!(pedcode,["P3/P4//P5/P6" for _ in 1:subpopsize])
    append!(nself,[4 for _ in 1:subpopsize])
    append!(sampleid,[string("RIL_",i) for i in 1:subpopsize])
    # 2ril-fixed: P5,P6
    append!(pedcode,["P5/P6" for _ in 1:subpopsize])
    append!(nself,["FIXED" for _ in 1:subpopsize])
    append!(sampleid,[string("FIXED_",i) for i in 1:subpopsize])
    peddf = DataFrame(sampleid=sampleid, pedcode=pedcode,nself=nself)
    CSV.write(getabsfile(workdir,outfile),peddf)
end
