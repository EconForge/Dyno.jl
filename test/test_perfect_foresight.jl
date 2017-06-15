
rootdir = Pkg.dir("Dyno")
filename = joinpath(rootdir,"examples","rbc.mod")

import Dyno
@time model = Dyno.import_model(filename, print_code=false)

model.functions.f_static

y = [model.calibration[p] for p in model.symbols[:endogenous]]
e = [model.calibration[p] for p in model.symbols[:exogenous]]
p = [model.calibration[p] for p in model.symbols[:parameters]]
using Dolang
model.functions.f_static(Der{1},[y;e],p)
model.functions.f_dynamic(Der{1},[y;y;y;e],p)



exo = [0.0, 0.1, 0.05, 0.0]''

sim = Dyno.deterministic(model, exo)


dr = Dyno.solve(model)
