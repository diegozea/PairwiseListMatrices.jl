isdefined(Base, :__precompile__) && __precompile__()

module PairwiseListMatrices

  using NamedArrays

  export @iterateupper,
  @iteratelist,
  @iteratediag,

  PairwiseListMatrix,
  hasdiagonal,
  getlist, getdiag,
  lengthlist, ij2k,
  sum_nodiag, mean_nodiag,
  zscore!, zscore,
  from_table, to_table

  # protovis

  "Function from MLPlots.jl, written by Tom Breloff."
  function _is_installed(name::String)
      try
          Pkg.installed(name) === nothing ? false: true
      catch
          false
      end
  end

  include("macros.jl")
  include("pairwiselistmatrix.jl")
  # include("protovis.jl")

  if _is_installed("DataFrames")
      include(joinpath("DataFrames","dataframe.jl"))
  end

end
