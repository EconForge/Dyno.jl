module Dyno

    export import_modfile, parse_modfile, import_data, solve, simulate, deterministic

    using Dolang
    using DataFrames
    using DataStructures
    using Distributions
    using Dolang
    import JSON
    import Dolang
    import Dolang: Der

    include("trisystem.jl")
    include("modfile_parsing.jl")
    include("model_import.jl")
    include("perturbation.jl")
    include("deterministic.jl")

end # module
