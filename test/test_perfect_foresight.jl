
rootdir = Pkg.dir("Dyno")
filename = joinpath(rootdir,"examples","rbc.mod")

model = Dyno.import_modfile(filename)

exo = [0.0, 0.1, 0.05, 0.0]''

sim = Dyno.deterministic(model, exo)
