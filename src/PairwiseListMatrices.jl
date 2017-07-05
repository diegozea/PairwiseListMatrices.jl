isdefined(Base, :__precompile__) && __precompile__()

module PairwiseListMatrices

"Function from MLPlots.jl, written by Tom Breloff."
function _is_installed(name::String)
    try
        Pkg.installed(name) === nothing ? false: true
    catch
        false
    end
end

if _is_installed("DataFrames")
    using DataFrames
    import DataFrames: join, zscore, zscore!
end

using NamedArrays
using RecipesBase

export  @iterateupper,
        @iteratelist,
        @iteratediag,

        PairwiseListMatrix,
        hasdiagonal,
        getlist, getdiag,
        lengthlist, ij2k,
        sum_nodiag, mean_nodiag,
        zscore!, zscore,
        getlabels,
        setlabels, setlabels!,
        from_table, to_table

include("macros.jl")
include("pairwiselistmatrix.jl")
include("plotrecipes.jl")

if _is_installed("DataFrames")
    include(joinpath("DataFrames","dataframe.jl"))
end

end
