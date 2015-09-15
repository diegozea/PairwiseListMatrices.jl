abstract AbstractPairwiseList{T, L} <: AbstractMatrix{T}
# AbstractPairwiseList fields: list, labels and nelements

abstract AbstractPairwiseListDiagonalMatrix{T, L} <: AbstractPairwiseList{T, L}
# Same fields than AbstractPairwiseList

abstract AbstractPairwiseListMatrix{T, L} <: AbstractPairwiseList{T, L}
# AbstractPairwiseListMatrix has also a diagonal field

type PairwiseListDiagonalSymmetric{T, L} <: AbstractPairwiseListDiagonalMatrix{T, L}
  list::Vector{T}
  labels::IndexedArray{L}
  nelements::Int
end

type PairwiseListSymmetric{T, L} <: AbstractPairwiseListMatrix{T, L}
  list::Vector{T}
  diagonal::Vector{T}
  labels::IndexedArray{L}
  nelements::Int
end

type PairwiseListDiagonalSquareTriangular{T, L} <: AbstractPairwiseListDiagonalMatrix{T, L}
  list::Vector{T}
  labels::IndexedArray{L}
  nelements::Int
  upper::Bool
end

type PairwiseListSquareTriangular{T, L} <: AbstractPairwiseListMatrix{T, L}
  list::Vector{T}
  diagonal::Vector{T}
  labels::IndexedArray{L}
  nelements::Int
  upper::Bool
end

# Creation
# ========

@inline _list_length(nelements) = @fastmath div(nelements*(nelements-1),2)
@inline _list_with_diagonal_length(nelements) = @fastmath div(nelements*(nelements+1),2)

_test_label_length(labels, nelements) = length(labels) == nelements || throw(ArgumentError("labels must have $nelements elements!"))

# Empties
# -------

function PairwiseListDiagonalSymmetric{T}(::Type{T}, nelements::Int, labels::AbstractVector=collect(1:nelements))
  _test_label_length(labels, nelements)
  PairwiseListDiagonalSymmetric(Array(T, _list_with_diagonal_length(nelements)), IndexedArray(labels), nelements)
end
PairwiseListDiagonalSymmetric{T}(::Type{T}, labels::AbstractVector) = PairwiseListDiagonalSymmetric(Array(T, _list_with_diagonal_length(length(labels))), IndexedArray(labels), length(labels))

function PairwiseListSymmetric{T}(::Type{T}, nelements::Int, labels::AbstractVector=collect(1:nelements))
  _test_label_length(labels, nelements)
  PairwiseListSymmetric(Array(T, _list_length(nelements)), Array(T, nelements), IndexedArray(labels), nelements)
end
PairwiseListSymmetric{T}(::Type{T}, labels::AbstractVector) = PairwiseListSymmetric(Array(T, _list_length(length(labels))), Array(T, length(labels)), IndexedArray(labels), length(labels))

function PairwiseListDiagonalSquareTriangular{T}(::Type{T}, nelements::Int, upper::Bool=true, labels::AbstractVector=collect(1:nelements))
  _test_label_length(labels, nelements)
  PairwiseListDiagonalSquareTriangular(Array(T, _list_with_diagonal_length(nelements)), IndexedArray(labels), nelements, upper)
end
function PairwiseListDiagonalSquareTriangular{T}(::Type{T}, upper::Bool, labels::AbstractVector)
  len = length(labels)
  PairwiseListDiagonalSquareTriangular(Array(T, _list_with_diagonal_length(len)), IndexedArray(labels), len, upper)
end

function PairwiseListSquareTriangular{T}(::Type{T}, nelements::Int, upper::Bool=true, labels::AbstractVector=collect(1:nelements))
  _test_label_length(labels, nelements)
  PairwiseListSquareTriangular(Array(T, _list_length(nelements)), Array(T, nelements), IndexedArray(labels), nelements, upper)
end
function PairwiseListSquareTriangular{T}(::Type{T}, upper::Bool, labels::AbstractVector)
  len = length(labels)
  PairwiseListSquareTriangular(Array(T, _list_length(len)), Array(T, len), IndexedArray(labels), len, upper)
end

# From a list
# -----------

@inline _nelements(len::Int) = @fastmath div(1+Int(sqrt(1+8*len)),2)
@inline _nelements_with_diagonal(len::Int) = @fastmath div(Int(sqrt(1+8*len)-1),2)

_test_diagonal_length(diagonal, nelements) = length(diagonal) == nelements || throw(ArgumentError("diagonal must have $nelements elements!"))

_set_diagonal{T <: Number}(value::T, nelements) = fill!(Array(T, nelements), value)
function _set_diagonal{T <: Number}(list::AbstractVector{T}, nelements)
  _test_diagonal_length(list, nelements)
  list
end

function PairwiseListDiagonalSymmetric(list::AbstractVector, labels::AbstractVector=collect(1:_nelements_with_diagonal(length(list))))
  nelements = _nelements_with_diagonal(length(list))
  _test_label_length(labels, nelements)
  PairwiseListDiagonalSymmetric(list, IndexedArray(labels), nelements)
end

function PairwiseListSymmetric{T}(list::AbstractVector{T}, labels::AbstractVector=collect(1:_nelements(length(list))), diagonal=zero(T))
  nelements = _nelements(length(list))
  _test_label_length(labels, nelements)
  PairwiseListSymmetric(list, _set_diagonal(diagonal, nelements), IndexedArray(labels), nelements)
end

function PairwiseListDiagonalSquareTriangular(list::AbstractVector, upper::Bool=true, labels::AbstractVector=collect(1:_nelements_with_diagonal(length(list))))
  nelements = _nelements_with_diagonal(length(list))
  _test_label_length(labels, nelements)
  PairwiseListDiagonalSquareTriangular(list, IndexedArray(labels), nelements, upper)
end

function PairwiseListSquareTriangular{T}(list::AbstractVector{T}, upper::Bool=true, labels::AbstractVector=collect(1:_nelements(length(list))), diagonal=zero(T))
  nelements = _nelements(length(list))
  _test_label_length(labels, nelements)
  PairwiseListSquareTriangular(list, _set_diagonal(diagonal, nelements), IndexedArray(labels), nelements, upper)
end

# Similar and Copy
# ================

for fun in (:similar, :copy)
  @eval begin
    $(fun){T,L}(lm::PairwiseListDiagonalSymmetric{T,L}) = PairwiseListDiagonalSymmetric{T,L}($(fun)(lm.list), copy(lm.labels), copy(lm.nelements))
    $(fun){T,L}(lm::PairwiseListSymmetric{T,L}) = PairwiseListSymmetric{T,L}($(fun)(lm.list), $(fun)(lm.diagonal), copy(lm.labels), copy(lm.nelements))
    $(fun){T,L}(lm::PairwiseListDiagonalSquareTriangular{T,L}) = PairwiseListDiagonalSquareTriangular{T,L}($(fun)(lm.list), copy(lm.labels), copy(lm.nelements), copy(lm.upper))
    $(fun){T,L}(lm::PairwiseListSquareTriangular{T,L}) = PairwiseListSquareTriangular{T,L}($(fun)(lm.list), $(fun)(lm.diagonal), copy(lm.labels), copy(lm.nelements), copy(lm.upper))
  end
end

similar{T,L,S}(lm::PairwiseListDiagonalSymmetric{T,L}, ::Type{S}) = PairwiseListDiagonalSymmetric{S,L}(similar(lm.list, S), copy(lm.labels), copy(lm.nelements))
similar{T,L,S}(lm::PairwiseListSymmetric{T,L}, ::Type{S}) = PairwiseListSymmetric{S,L}(similar(lm.list, S), similar(lm.diagonal, S), copy(lm.labels), copy(lm.nelements))
similar{T,L,S}(lm::PairwiseListDiagonalSquareTriangular{T,L}, ::Type{S}) = PairwiseListDiagonalSquareTriangular{S,L}(similar(lm.list, S), copy(lm.labels), copy(lm.nelements), copy(lm.upper))
similar{T,L,S}(lm::PairwiseListSquareTriangular{T,L}, ::Type{S}) = PairwiseListSquareTriangular{S,L}(similar(lm.list, S), similar(lm.diagonal, S), copy(lm.labels), copy(lm.nelements), copy(lm.upper))

# Common methods
# ==============

size(lm::AbstractPairwiseList) = (lm.nelements, lm.nelements)
length(lm::AbstractPairwiseList) = lm.nelements * lm.nelements

eltype{T, L}(lm::AbstractPairwiseList{T, L}) = T

# Indexing (getindex)
# ===================

@inline _listindex(i, j, n) = @fastmath div((n*(n-1))-((n-i)*(n-i-1)),2) - n + j
@inline _listindex_with_diagonal(i, j, n) = @fastmath div((n*(n+1))-((n-i)*(n-i+1)),2) - n + j

function getindex(lm::PairwiseListSymmetric, i::Int, j::Int)
  if i < j
    return(lm.list[_listindex(i, j, lm.nelements)])
  elseif i > j
    return(lm.list[_listindex(j, i, lm.nelements)])
  else
    return(lm.diagonal[i])
  end
end

getindex(lm::PairwiseListDiagonalSymmetric, i::Int, j::Int) = i <= j ? lm.list[_listindex_with_diagonal(i, j, lm.nelements)] : lm.list[_listindex_with_diagonal(j, i, lm.nelements)]

function getindex{T, L}(lm::PairwiseListSquareTriangular{T, L}, i::Int, j::Int)
  if i == j
    return(lm.diagonal[i])
  end
  if lm.upper && i < j
    return(lm.list[_listindex(i, j, lm.nelements)])
  elseif !lm.upper && i > j
    return(lm.list[_listindex(j, i, lm.nelements)])
  else
    return(zero(T))
  end
end

function getindex{T, L}(lm::PairwiseListDiagonalSquareTriangular{T, L}, i::Int, j::Int)
  if lm.upper && i <= j
    return(lm.list[_listindex_with_diagonal(i, j, lm.nelements)])
  elseif !lm.upper && i >= j
    return(lm.list[_listindex_with_diagonal(j, i, lm.nelements)])
  else
    return(zero(T))
  end
end

# Indexing using labels (getlabel)
# ================================

getlabel(lm::AbstractPairwiseList, i, j) =  getindex(lm, findfirst(lm.labels, i), findfirst(lm.labels, j))

# Set values (setindex!)
# ======================

function setindex!(lm::PairwiseListSymmetric, v, i::Int, j::Int)
  if i < j
    return(setindex!(lm.list, v, _listindex(i, j, lm.nelements)))
  elseif i > j
    return(setindex!(lm.list, v, _listindex(j, i, lm.nelements)))
  else
    return(setindex!(lm.diagonal, v, i))
  end
end

setindex!(lm::PairwiseListDiagonalSymmetric, v, i::Int, j::Int) = i <= j ? setindex!(lm.list, v, _listindex_with_diagonal(i, j, lm.nelements)) :
  setindex!(lm.list, v, _listindex_with_diagonal(j, i, lm.nelements))

function setindex!(lm::PairwiseListSquareTriangular, v, i::Int, j::Int)
  if i == j
    return(setindex!(lm.diagonal, v, i))
  end
  if lm.upper && i < j
    return(setindex!(lm.list, v, _listindex(i, j, lm.nelements)))
  elseif !lm.upper && i > j
    return(setindex!(lm.list, v, _listindex(j, i, lm.nelements)))
  else
    throw(BoundsError(lm,(i,j)))
  end
end

function setindex!(lm::PairwiseListDiagonalSquareTriangular, v, i::Int, j::Int)
  if lm.upper && i <= j
    return(setindex!(lm.list, v, _listindex_with_diagonal(i, j, lm.nelements)))
  elseif !lm.upper && i >= j
    return(setindex!(lm.list, v, _listindex_with_diagonal(j, i, lm.nelements)))
  else
    throw(BoundsError(lm,(i,j)))
  end
end

# Set values using labels (setlabel!)
# ===================================

setlabel!(lm::AbstractPairwiseList, v, i, j) =  setindex!(lm, v, findfirst(lm.labels, i), findfirst(lm.labels, j))

# Transpose
# =========

transpose(lm::Union(PairwiseListDiagonalSymmetric, PairwiseListSymmetric)) = lm
transpose!(lm::Union(PairwiseListDiagonalSymmetric, PairwiseListSymmetric)) = lm
ctranspose(lm::Union(PairwiseListDiagonalSymmetric, PairwiseListSymmetric)) = lm
ctranspose!(lm::Union(PairwiseListDiagonalSymmetric, PairwiseListSymmetric)) = lm

transpose!(lm::Union(PairwiseListDiagonalSquareTriangular, PairwiseListSquareTriangular)) = (lm.upper = !lm.upper ; lm)
transpose(lm::Union(PairwiseListDiagonalSquareTriangular, PairwiseListSquareTriangular)) = (mat = copy(lm); mat.upper = !mat.upper ; mat)
ctranspose!(lm::Union(PairwiseListDiagonalSquareTriangular, PairwiseListSquareTriangular)) = (lm.upper = !lm.upper ; lm)
ctranspose(lm::Union(PairwiseListDiagonalSquareTriangular, PairwiseListSquareTriangular)) = (mat = copy(lm); mat.upper = !mat.upper ; mat)

# diag and full
# =============

diag(lm::AbstractPairwiseListMatrix) = lm.diagonal

full{T, L}(lm::AbstractPairwiseList{T, L}) = T[ lm[i,j]  for i in 1:lm.nelements, j in 1:lm.nelements ]


# Unary operations
# ================

for una in (:abs, :-)
  @eval begin
    $(una){T,L}(lm::PairwiseListDiagonalSymmetric{T,L}) = PairwiseListDiagonalSymmetric{T,L}($(una)(lm.list), lm.labels, lm.nelements)
    $(una){T,L}(lm::PairwiseListSymmetric{T,L}) = PairwiseListSymmetric{T,L}($(una)(lm.list), $(una)(lm.diagonal), lm.labels, lm.nelements)
    $(una){T,L}(lm::PairwiseListDiagonalSquareTriangular{T,L}) = PairwiseListDiagonalSquareTriangular{T,L}($(una)(lm.list), lm.labels, lm.nelements, lm.upper)
    $(una){T,L}(lm::PairwiseListSquareTriangular{T,L}) = PairwiseListSquareTriangular{T,L}($(una)(lm.list), $(una)(lm.diagonal), lm.labels, lm.nelements, lm.upper)
  end
end

# Binary operations
# =================

for bin in ( :-, :+, :.*, :./, :.+, :.- )

  @eval begin

    function $(bin){T,L}(A::PairwiseListDiagonalSymmetric{T,L}, B::PairwiseListDiagonalSymmetric{T,L})
      if A.labels != B.labels || A.nelements != B.nelements
        return($(bin)(full(A), full(B)))
      end
      PairwiseListDiagonalSymmetric{T,L}($(bin)(A.list, B.list), A.labels, A.nelements)
    end

    function $(bin){T,L}(A::PairwiseListSymmetric{T,L}, B::PairwiseListSymmetric{T,L})
      if A.labels != B.labels || A.nelements != B.nelements
        return($(bin)(full(A), full(B)))
      end
      PairwiseListSymmetric{T,L}($(bin)(A.list, B.list), $(bin)(A.diagonal, B.diagonal), A.labels, A.nelements)
    end

    function $(bin){T,L}(A::PairwiseListDiagonalSquareTriangular{T,L}, B::PairwiseListDiagonalSquareTriangular{T,L})
      if A.labels != B.labels || A.nelements != B.nelements || A.upper != B.upper
        return($(bin)(full(A), full(B)))
      end
      PairwiseListDiagonalSquareTriangular{T,L}($(bin)(A.list, B.list), A.labels, A.nelements, A.upper)
    end

    function $(bin){T,L}(A::PairwiseListSquareTriangular{T,L}, B::PairwiseListSquareTriangular{T,L})
      if A.labels != B.labels || A.nelements != B.nelements || A.upper != B.upper
        return($(bin)(full(A), full(B)))
      end
      PairwiseListSquareTriangular{T,L}($(bin)(A.list, B.list), $(bin)(A.diagonal, B.diagonal), A.labels, A.nelements, A.upper)
    end

  end

end
