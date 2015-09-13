module PairedListMatrices

  using IndexedArrays

  import Base: size, getindex, setindex!, length, eltype, similar, copy, abs, full,
  transpose, transpose!, ctranspose, ctranspose!, diag, -, +, .*, .-, .+, .-, ./

  export AbstractPairedList, AbstractPairedListDiagonalMatrix, AbstractPairedListMatrix,
  PairedListDiagonalSymmetric, PairedListSymmetric, PairedListDiagonalSquareTriangular, PairedListSquareTriangular,
  getlabel, setlabel!

  include("matrices.jl")

end # module
