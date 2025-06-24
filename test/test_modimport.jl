# using CondaPkg

# preproc_exe = joinpath(CondaPkg.bindirs()[1], "dynare-preprocessor")

# json_dump = 

import Dyno
using StaticArrays

s1 = time()


function compute()

    model = Dyno.import_modfile("examples/example1.mod")

    p = SVector(model.calibration_grouped[:parameters]...)
    e = SVector(model.calibration_grouped[:exogenous]...)
    y = SVector(model.calibration_grouped[:endogenous]...)

    r = Dyno.static(model, y, e, p)
    return r

end



model = Dyno.import_modfile("examples/example1.mod")

p = SVector(model.calibration_grouped[:parameters]...)
e = SVector(model.calibration_grouped[:exogenous]...)
y = SVector(model.calibration_grouped[:endogenous]...)

r = Dyno.static(model, y, e, p)