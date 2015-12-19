isdefined(Base, :__precompile__) && __precompile__()

module PairwiseListMatrices

  using IndexedArrays
  using Requires

  import Base: size, getindex, setindex!, length, eltype,
  similar, copy, abs, full, zeros, ones, convert,
  transpose, transpose!, ctranspose, ctranspose!, diag, svd,
  triu!, triu, writedlm, writecsv,
  -, +, .*, .-, .+, .-, ./, *, /,
  # faster
  mean, sum, varm, var, std, sqrt

  export @iterateupper,
  @iteratelist,
  @iteratediag,

  PairwiseListMatrix, hasdiagonal, lengthlist,
  getlist, getdiag, labels, labels!, getlabel, setlabel!,
  sum_nodiag, mean_nodiag, zscore!, zscore,
  from_table, to_table,

  protovis

  include("macros.jl")
  include("pairwiselistmatrix.jl")
  include("protovis.jl")

  @require DataFrames include(joinpath("DataFrames","dataframe.jl"))

end
