module PairwiseListMatrices

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
        to_dataframe, from_dataframe

include("macros.jl")
include("pairwiselistmatrix.jl")
include("plotrecipes.jl")

function to_dataframe(args...)
    throw(ErrorException("Deprecated function, use DataFrame(to_dict(args...)) instead."))
end

function from_dataframe(args...)
    throw(ErrorException("Deprecated function, use from_table(args...) instead."))
end

end
