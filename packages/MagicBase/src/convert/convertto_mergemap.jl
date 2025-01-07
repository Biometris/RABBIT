function tomergemapfile(mapfile::AbstractString;
    chrsubset::Union{Nothing,AbstractVector} = nothing, 
    missingstring=["NA","missing"],    
    commentstring::AbstractString="##",
    outstem::AbstractString="outstem",
    workdir::AbstractString = pwd())
    inputmap = readmarkermap(mapfile; missingstring, del_ungrouped = true, commentstring, workdir)
    ginputmap = groupby(inputmap,:linkagegroup)    
    if !isnothing(chrsubset) && eltype(chrsubset) <: Integer
         nchr  = length(ginputmap)
         b = 1 .<= chrsubset .<= nchr
         ginputmap = ginputmap[chrsubset[b]]
    end
    outfile = getabsfile(workdir, outstem*".txt")
    open(outfile,"w") do io
        for df in ginputmap
            lgname = df[1,:linkagegroup]
            write(io,"group\t",lgname,"\n")
            write(io, ";BEGINOFGROUP\n")
            for row in eachrow(df)
                write(io, string(row[:marker],"\t",row[:poscm]),"\n")
            end
            write(io, ";ENDOFGROUP\n")
        end
    end
    outfile
end