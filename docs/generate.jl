execute = isempty(ARGS) || ARGS[1] == "run"

using Literate

lit = joinpath(@__DIR__, "lit")
src = joinpath(@__DIR__, "src")
gen = joinpath(@__DIR__, "src/generated")

for (root, _, files) in walkdir(lit), file in files
    splitext(file)[2] == ".jl" || continue # process .jl files only
    ipath = joinpath(root, file)
    opath = splitdir(replace(ipath, lit => gen))[1]
    Literate.markdown(ipath, opath, documenter = execute)
    Literate.notebook(ipath, opath; execute = false)
end

# functions
ismd(f) = splitext(f)[2] == ".md"
pages(folder) =
    [joinpath("generated/", folder, f) for f in readdir(joinpath(gen, folder)) if ismd(f)]
