

function plotmarkererror(mapfile::AbstractString;
    tukeyfence::Real=3.0,
    missingstring=["NA","missing"],
    commentstring::AbstractString="##",
    workdir::AbstractString = pwd())
    markermap = readmarkermap(mapfile; del_ungrouped = false,        
        commentstring,missingstring, workdir)  
    plotmarkererror(markermap;tukeyfence)
end

function plotmarkererror(markermap::AbstractDataFrame;
    tukeyfence::Real=3.0)
    colls = [:foundererror,:offspringerror,:seqerror,:allelebalancemean,:allelebalancedisperse,:alleledropout]
    fencels = tukeyfence*ones(length(colls))
    b = [begin         
        isanymiss = any(ismissing.(markermap[!,i]))         
        if isanymiss 
            true
        else
            ls = unique(round.(unique(markermap[!,i]),digits=4))
            # length(ls) <=1 && println("col=",i, ",unqiue values=",ls)
            length(ls) <=1 
        end
    end for i in colls]
    deleteat!(colls,b)
    deleteat!(fencels,b)
    isempty(colls) && return nothing    
    gls = [begin     
        col = colls[c]        
        whisker_range = fencels[c]        
        # row = [!(ismissing(i) || isnan(i)) for i in markermap[!,col]]
        row = !
        yls = Vector(markermap[row, col])
        yls[yls .== 0.0] .= 1e-6
        if in(col, [:foundererror,:offspringerror,:seqerror,:allelebalancedisperse,:alleledropout]) 
            @. yls = log10(yls)            
            ylabel = string("log10 ", col)
        else
            # @. yls = log10(abs(yls))
            # ylabel = string("log10 |", col,"|")            
            ylabel = string(col)
        end                     
        xls = string.(markermap[row,:linkagegroup])                
        g =violin(xls, yls;
            fillalpha=0.75, linewidth=0.5,                
            legend=false,            
            ylabel
        )              
        boxplot!(xls, yls;
            fillalpha=0.75, linewidth=2,
            whisker_range, 
            legend=false, 
            title = string(col, ": tukeyfence=",whisker_range)
        )   
        g
    end for c in eachindex(colls)]
    nchr = length(unique(markermap[!,:linkagegroup]))
    plot(gls...;
        layout=(length(gls),1),
        size = (max(1000,nchr*60),300*length(gls)),
        left_margin=20Plots.mm,         
        bottom_margin=15Plots.mm,
    )
end


function plot_peroffspringerror(perofffile::AbstractString;
    tukeyfence::Real=3.0,    
    workdir::AbstractString = pwd())
    peroffspringerr = CSV.read(getabsfile(workdir, perofffile),DataFrame)    
    plot_peroffspringerror(peroffspringerr; tukeyfence)
end


function plot_peroffspringerror(peroffspringerr::AbstractDataFrame;
    tukeyfence::Real=3.0)
    whisker_range = tukeyfence
    colnames = names(peroffspringerr)
    colls = findall(occursin.(r"^peroffspringerr",colnames))
    yls = vec(Matrix(peroffspringerr[!,colls]))
    yls[yls .< 1e-6] .= 1e-6
    yls .= log10.(yls)
    labs = replace.(colnames[colls],"peroffspringerror_"=>"")
    xls = repeat(labs,inner=size(peroffspringerr,1))
    fig =violin(xls, yls;
        fillalpha=0.75, linewidth=0.5,                
        legend=false,            
        ylabel = "log10 peroffspringerr"
    )              
    boxplot!(xls, yls;
        fillalpha=0.75, linewidth=2,            
        legend=false, 
        title = string("tukeyfence=",whisker_range)
    )   
    fig
end