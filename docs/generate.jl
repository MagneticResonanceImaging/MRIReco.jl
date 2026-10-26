execute = isempty(ARGS) || ARGS[1] == "run"

using Literate

# Literate sources in src/literate become Markdown pages and notebooks in src/generated,
# with the same subfolders (tutorials, howto)
lit = joinpath(@__DIR__, "src/literate")
gen = joinpath(@__DIR__, "src/generated")

for (root, _, files) in walkdir(lit), file in files
    splitext(file)[2] == ".jl" || continue # process .jl files only
    ipath = joinpath(root, file)
    opath = splitdir(replace(ipath, lit => gen))[1]
    Literate.markdown(ipath, opath, documenter = execute)
    Literate.notebook(ipath, opath; execute = false)
end
