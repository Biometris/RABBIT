module HMM

using Distributions: Categorical

include("posteriordecode.jl")

include("posteriorsample.jl")

include("logforback.jl")

include("viterbidecode.jl")

include("pathformat.jl")

end # module
