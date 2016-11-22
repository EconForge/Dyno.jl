
rootdir = Pkg.dir("Dyno")
filename = joinpath(rootdir,"examples","rbc.mod")

import Dolang


import Dyno
model = Dyno.import_modfile(filename)

sol = Dyno.solve(model)
Dyno.simulate(sol)
