# Dyno

(very preliminary version)

## Install

Dyno requires julia 0.6 and depends on the following Julia packages:

- Dolang: `Pkg.add("Dolang")`
- DataFrames: `Pkg.add("DataFrames")`

- Plots (optional): `Pkg.add("Plots") (or any other plotting library)

Then:
- `Pkg.clone("https://github.com/EconForge/Dyno.jl.git")`

Dyno depends on the dynare executable. It is downloaded (linux only) with:
- `Pkg.build("Dyno")`


## Example

From `example.jl`:

```
#####
##### stochastic simulation
#####

using Dyno

rootdir = Pkg.dir("Dyno")
filename = joinpath(rootdir,"examples","rbc.mod")

@time model = Model(filename)
@time sol = solve(model)
@time sim = simulate(sol)
@time sim = simulate(model)#

using Plots
plot(sim[:t], sim[:y])

#####
##### deterministic solve
#####

exo  = [ 0.0 0.1 0.00]'
@time sol = deterministic(model, exo, verbose=true,N=200)

using Plots
plot(sol[:t], sol[:y])
```
