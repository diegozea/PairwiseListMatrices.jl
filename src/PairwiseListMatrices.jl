__precompile__()

"""
This package allows you to use a pairwise list as a symmetric matrix.

```julia
using PairwiseListMatrices
```
"""
module PairwiseListMatrices

# standard libraries
using DelimitedFiles
using LinearAlgebra
using Statistics

# external packages
using NamedArrays
using RecipesBase

export  @iterateupper,
        @iteratelist,
        @iteratediag,

        PairwiseListMatrix,
        hasdiagonal,
        getlist, getdiag,
        diagonal,
        lengthlist, ij2k,
        sum_nodiag, mean_nodiag,
        zscore!, zscore,
        getlabels,
        setlabels, setlabels!,
        from_table, to_table,
        to_dict,
        join,
        writedlm

include("macros.jl")
include("pairwiselistmatrix.jl")
include("plotrecipes.jl")

end
