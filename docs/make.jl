using Documenter, PairwiseListMatrices

makedocs(
    doctest = true,
    format = :html,
    sitename = "PairwiseListMatrices",
    modules = [PairwiseListMatrices],
    pages = [
        "index.md",
        "api.md"
    ]
)

deploydocs(
    repo   = "github.com/diegozea/PairwiseListMatrices.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing
)
