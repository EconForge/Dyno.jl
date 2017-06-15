# Dyno

(very preliminary version)

## Install

Dyno requires julia 0.5 and depends on the following Julia packages:

- Dolang: `Pkg.add("Dolang")`
- DataFrames: `Pkg.add("DataFrames")`

- SymEngine (optional): `Pkg.add("SymEngine") (for faster calculations)
- Gadfly (optional): `Pkg.add("Gadfly") (or any other plotting library)

## Example

From `example.jl`:

```
#####
##### stochastic simulation
#####

using Dyno

rootdir = Pkg.dir("Dyno")
filename = joinpath(rootdir,"examples","rbc.mod")

@time model = import_modfile(filename)
@time sol = solve(model)
@time sim = simulate(sol)
@time sim = simulate(model)#

using Gadfly
plot(x=sim[:t], y=sim[:y], Geom.line)

#####
##### deterministic solve
#####

exo  = [ 0.0 0.1 0.00]'
@time sol = deterministic(model, exo, verbose=true,N=200)

using Gadfly
plot(x=sol[:t], y=sol[:y], Geom.line)
```
