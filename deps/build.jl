if is_linux()
    path = joinpath(Pkg.dir("Dyno", "deps", "usr", "bin"))
    mkpath(path)
    filename = joinpath(path, "dynare_m")
    download("https://github.com/albop/demo_dyno/raw/master/dynare_m", filename)
    run(`chmod +x $filename`)
end

# using BinDeps
# @BinDeps.setup
# dynare = library_dependency("dynare")
# provides(Binaries, URI("https://github.com/albop/demo_dyno/raw/master/dynare_m"), dynare, os=:linux)
# @BinDeps.install
