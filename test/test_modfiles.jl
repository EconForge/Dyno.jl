import Dyno

examples_dir = joinpath( Pkg.dir("Dyno"), "examples")


filename = joinpath(examples_dir, "rbc.mod")
@time model = Dyno.import_model(filename)
@time dr = Dyno.solve(model)

filename =  joinpath(examples_dir,"fs2000_working.mod")
@time model = Dyno.import_model(filename)
@time dr = Dyno.solve(model)

filename =  joinpath(examples_dir,"example1.mod")
@time model = Dyno.import_model(filename)
@time dr = Dyno.solve(model)


# NK_baseline currently needs a steady-state file (experimental)
include( joinpath(examples_dir,"NK_baseline_steadystate.jl") )
import NK_baseline

filename = joinpath(examples_dir, "NK_baseline.mod")
@time model = Dyno.import_model(filename, print_code=false)
new_calib = NK_baseline.steadystate(model)
NK_baseline.update_calibration(model, new_calib)
@time dr = Dyno.solve(model)
