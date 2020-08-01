using Documenter, PairwiseListMatrices

makedocs(
    doctest=true,
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", nothing) == "true"
    ),
    sitename="PairwiseListMatrices",
    modules=[PairwiseListMatrices],
    pages=[
        "index.md",
        "api.md"
    ]
)

deploydocs(
    repo="github.com/diegozea/PairwiseListMatrices.jl.git",
    target="build",
    deps=nothing,
    make=nothing
)
