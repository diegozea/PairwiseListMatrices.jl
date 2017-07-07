using Documenter, PairwiseListMatrices

makedocs(
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
    julia  = 0.5,
    deps   = nothing,
    make   = nothing
)
