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

# does not work
filename = joinpath(examples_dir, "EA_QUEST3_working.mod")
@time model = Dyno.import_model(filename)
@time dr = Dyno.solve(model)
