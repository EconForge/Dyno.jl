#####
##### stochastic simulation
#####

using Dyno

rootdir = Pkg.dir("Dyno")
filename = joinpath(rootdir,"examples","rbc.mod")

@time model = import_modfile(filename)
@time sol = solve(model)
@time sim = simulate(sol)
@time sim = simulate(model)

using Gadfly
plot(x=sim[:t], y=sim[:y], Geom.line)

#####
##### deterministic solve
#####

exo  = [ 0.0 0.1 0.00]'
@time sol = deterministic(model, exo, verbose=true,N=200)

using Gadfly
plot(x=sol[:t], y=sol[:y], Geom.line)
