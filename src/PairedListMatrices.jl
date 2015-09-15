module PairwiseListMatrices

  using IndexedArrays

  import Base: size, getindex, setindex!, length, eltype, similar, copy, abs, full,
  transpose, transpose!, ctranspose, ctranspose!, diag, -, +, .*, .-, .+, .-, ./

  export AbstractPairwiseList, AbstractPairwiseListDiagonalMatrix, AbstractPairwiseListMatrix,
  PairwiseListDiagonalSymmetric, PairwiseListSymmetric, PairwiseListDiagonalSquareTriangular, PairwiseListSquareTriangular,
  getlabel, setlabel!

  include("matrices.jl")

end # module
