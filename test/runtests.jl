using Dyno
using Base.Test

rootdir = Pkg.dir("Dyno")

include("test_model_import.jl")

include("test_trisystem.jl")
