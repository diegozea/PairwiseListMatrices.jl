abstract AbstractPairedList{T, L} <: AbstractMatrix{T}
# AbstractPairedList fields: list, labels and nelements

abstract AbstractPairedListDiagonalMatrix{T, L} <: AbstractPairedList{T, L}
# Same fields than AbstractPairedList

abstract AbstractPairedListMatrix{T, L} <: AbstractPairedList{T, L}
# AbstractPairedListMatrix has also a diagonal field

type PairedListDiagonalSymmetric{T, L} <: AbstractPairedListDiagonalMatrix{T, L}
  list::Vector{T}
  labels::IndexedArray{L}
  nelements::Int
end

type PairedListSymmetric{T, L} <: AbstractPairedListMatrix{T, L}
  list::Vector{T}
  diagonal::Vector{T}
  labels::IndexedArray{L}
  nelements::Int
end

type PairedListDiagonalSquareTriangular{T, L} <: AbstractPairedListDiagonalMatrix{T, L}
  list::Vector{T}
  labels::IndexedArray{L}
  nelements::Int
  upper::Bool
end

type PairedListSquareTriangular{T, L} <: AbstractPairedListMatrix{T, L}
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

function PairedListDiagonalSymmetric{T}(::Type{T}, nelements::Int, labels::AbstractVector=collect(1:nelements))
  _test_label_length(labels, nelements)
  PairedListDiagonalSymmetric(Array(T, _list_with_diagonal_length(nelements)), IndexedArray(labels), nelements)
end
PairedListDiagonalSymmetric{T}(::Type{T}, labels::AbstractVector) = PairedListDiagonalSymmetric(Array(T, _list_with_diagonal_length(length(labels))), IndexedArray(labels), length(labels))

function PairedListSymmetric{T}(::Type{T}, nelements::Int, labels::AbstractVector=collect(1:nelements))
  _test_label_length(labels, nelements)
  PairedListSymmetric(Array(T, _list_length(nelements)), Array(T, nelements), IndexedArray(labels), nelements)
end
PairedListSymmetric{T}(::Type{T}, labels::AbstractVector) = PairedListSymmetric(Array(T, _list_length(length(labels))), Array(T, length(labels)), IndexedArray(labels), length(labels))

function PairedListDiagonalSquareTriangular{T}(::Type{T}, nelements::Int, upper::Bool=true, labels::AbstractVector=collect(1:nelements))
  _test_label_length(labels, nelements)
  PairedListDiagonalSquareTriangular(Array(T, _list_with_diagonal_length(nelements)), IndexedArray(labels), nelements, upper)
end
function PairedListDiagonalSquareTriangular{T}(::Type{T}, upper::Bool, labels::AbstractVector)
  len = length(labels)
  PairedListDiagonalSquareTriangular(Array(T, _list_with_diagonal_length(len)), IndexedArray(labels), len, upper)
end

function PairedListSquareTriangular{T}(::Type{T}, nelements::Int, upper::Bool=true, labels::AbstractVector=collect(1:nelements))
  _test_label_length(labels, nelements)
  PairedListSquareTriangular(Array(T, _list_length(nelements)), Array(T, nelements), IndexedArray(labels), nelements, upper)
end
function PairedListSquareTriangular{T}(::Type{T}, upper::Bool, labels::AbstractVector)
  len = length(labels)
  PairedListSquareTriangular(Array(T, _list_length(len)), Array(T, len), IndexedArray(labels), len, upper)
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

function PairedListDiagonalSymmetric(list::AbstractVector, labels::AbstractVector=collect(1:_nelements_with_diagonal(length(list))))
  nelements = _nelements_with_diagonal(length(list))
  _test_label_length(labels, nelements)
  PairedListDiagonalSymmetric(list, IndexedArray(labels), nelements)
end

function PairedListSymmetric{T}(list::AbstractVector{T}, labels::AbstractVector=collect(1:_nelements(length(list))), diagonal=zero(T))
  nelements = _nelements(length(list))
  _test_label_length(labels, nelements)
  PairedListSymmetric(list, _set_diagonal(diagonal, nelements), IndexedArray(labels), nelements)
end

function PairedListDiagonalSquareTriangular(list::AbstractVector, upper::Bool=true, labels::AbstractVector=collect(1:_nelements_with_diagonal(length(list))))
  nelements = _nelements_with_diagonal(length(list))
  _test_label_length(labels, nelements)
  PairedListDiagonalSquareTriangular(list, IndexedArray(labels), nelements, upper)
end

function PairedListSquareTriangular{T}(list::AbstractVector{T}, upper::Bool=true, labels::AbstractVector=collect(1:_nelements(length(list))), diagonal=zero(T))
  nelements = _nelements(length(list))
  _test_label_length(labels, nelements)
  PairedListSquareTriangular(list, _set_diagonal(diagonal, nelements), IndexedArray(labels), nelements, upper)
end

# Similar and Copy
# ================

for fun in (:similar, :copy)
  @eval begin
    $(fun){T,L}(lm::PairedListDiagonalSymmetric{T,L}) = PairedListDiagonalSymmetric{T,L}($(fun)(lm.list), copy(lm.labels), copy(lm.nelements))
    $(fun){T,L}(lm::PairedListSymmetric{T,L}) = PairedListSymmetric{T,L}($(fun)(lm.list), $(fun)(lm.diagonal), copy(lm.labels), copy(lm.nelements))
    $(fun){T,L}(lm::PairedListDiagonalSquareTriangular{T,L}) = PairedListDiagonalSquareTriangular{T,L}($(fun)(lm.list), copy(lm.labels), copy(lm.nelements), copy(lm.upper))
    $(fun){T,L}(lm::PairedListSquareTriangular{T,L}) = PairedListSquareTriangular{T,L}($(fun)(lm.list), $(fun)(lm.diagonal), copy(lm.labels), copy(lm.nelements), copy(lm.upper))
  end
end

similar{T,L,S}(lm::PairedListDiagonalSymmetric{T,L}, ::Type{S}) = PairedListDiagonalSymmetric{S,L}(similar(lm.list, S), copy(lm.labels), copy(lm.nelements))
similar{T,L,S}(lm::PairedListSymmetric{T,L}, ::Type{S}) = PairedListSymmetric{S,L}(similar(lm.list, S), similar(lm.diagonal, S), copy(lm.labels), copy(lm.nelements))
similar{T,L,S}(lm::PairedListDiagonalSquareTriangular{T,L}, ::Type{S}) = PairedListDiagonalSquareTriangular{S,L}(similar(lm.list, S), copy(lm.labels), copy(lm.nelements), copy(lm.upper))
similar{T,L,S}(lm::PairedListSquareTriangular{T,L}, ::Type{S}) = PairedListSquareTriangular{S,L}(similar(lm.list, S), similar(lm.diagonal, S), copy(lm.labels), copy(lm.nelements), copy(lm.upper))

# Common methods
# ==============

size(lm::AbstractPairedList) = (lm.nelements, lm.nelements)
length(lm::AbstractPairedList) = lm.nelements * lm.nelements

eltype{T, L}(lm::AbstractPairedList{T, L}) = T

# Indexing (getindex)
# ===================

@inline _listindex(i, j, n) = @fastmath div((n*(n-1))-((n-i)*(n-i-1)),2) - n + j
@inline _listindex_with_diagonal(i, j, n) = @fastmath div((n*(n+1))-((n-i)*(n-i+1)),2) - n + j

function getindex(lm::PairedListSymmetric, i::Int, j::Int)
  if i < j
    return(lm.list[_listindex(i, j, lm.nelements)])
  elseif i > j
    return(lm.list[_listindex(j, i, lm.nelements)])
  else
    return(lm.diagonal[i])
  end
end

getindex(lm::PairedListDiagonalSymmetric, i::Int, j::Int) = i <= j ? lm.list[_listindex_with_diagonal(i, j, lm.nelements)] : lm.list[_listindex_with_diagonal(j, i, lm.nelements)]

function getindex{T, L}(lm::PairedListSquareTriangular{T, L}, i::Int, j::Int)
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

function getindex{T, L}(lm::PairedListDiagonalSquareTriangular{T, L}, i::Int, j::Int)
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

getlabel(lm::AbstractPairedList, i, j) =  getindex(lm, findfirst(lm.labels, i), findfirst(lm.labels, j))

# Set values (setindex!)
# ======================

function setindex!(lm::PairedListSymmetric, v, i::Int, j::Int)
  if i < j
    return(setindex!(lm.list, v, _listindex(i, j, lm.nelements)))
  elseif i > j
    return(setindex!(lm.list, v, _listindex(j, i, lm.nelements)))
  else
    return(setindex!(lm.diagonal, v, i))
  end
end

setindex!(lm::PairedListDiagonalSymmetric, v, i::Int, j::Int) = i <= j ? setindex!(lm.list, v, _listindex_with_diagonal(i, j, lm.nelements)) :
  setindex!(lm.list, v, _listindex_with_diagonal(j, i, lm.nelements))

function setindex!(lm::PairedListSquareTriangular, v, i::Int, j::Int)
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

function setindex!(lm::PairedListDiagonalSquareTriangular, v, i::Int, j::Int)
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

setlabel!(lm::AbstractPairedList, v, i, j) =  setindex!(lm, v, findfirst(lm.labels, i), findfirst(lm.labels, j))

# Transpose
# =========

transpose(lm::Union(PairedListDiagonalSymmetric, PairedListSymmetric)) = lm
transpose!(lm::Union(PairedListDiagonalSymmetric, PairedListSymmetric)) = lm
ctranspose(lm::Union(PairedListDiagonalSymmetric, PairedListSymmetric)) = lm
ctranspose!(lm::Union(PairedListDiagonalSymmetric, PairedListSymmetric)) = lm

transpose!(lm::Union(PairedListDiagonalSquareTriangular, PairedListSquareTriangular)) = (lm.upper = !lm.upper ; lm)
transpose(lm::Union(PairedListDiagonalSquareTriangular, PairedListSquareTriangular)) = (mat = copy(lm); mat.upper = !mat.upper ; mat)
ctranspose!(lm::Union(PairedListDiagonalSquareTriangular, PairedListSquareTriangular)) = (lm.upper = !lm.upper ; lm)
ctranspose(lm::Union(PairedListDiagonalSquareTriangular, PairedListSquareTriangular)) = (mat = copy(lm); mat.upper = !mat.upper ; mat)

# diag and full
# =============

diag(lm::AbstractPairedListMatrix) = lm.diagonal

full{T, L}(lm::AbstractPairedList{T, L}) = T[ lm[i,j]  for i in 1:lm.nelements, j in 1:lm.nelements ]


# Unary operations
# ================

for una in (:abs, :-)
  @eval begin
    $(una){T,L}(lm::PairedListDiagonalSymmetric{T,L}) = PairedListDiagonalSymmetric{T,L}($(una)(lm.list), lm.labels, lm.nelements)
    $(una){T,L}(lm::PairedListSymmetric{T,L}) = PairedListSymmetric{T,L}($(una)(lm.list), $(una)(lm.diagonal), lm.labels, lm.nelements)
    $(una){T,L}(lm::PairedListDiagonalSquareTriangular{T,L}) = PairedListDiagonalSquareTriangular{T,L}($(una)(lm.list), lm.labels, lm.nelements, lm.upper)
    $(una){T,L}(lm::PairedListSquareTriangular{T,L}) = PairedListSquareTriangular{T,L}($(una)(lm.list), $(una)(lm.diagonal), lm.labels, lm.nelements, lm.upper)
  end
end

# Binary operations
# =================

for bin in ( :-, :+, :.*, :./, :.+, :.- )

  @eval begin

    function $(bin){T,L}(A::PairedListDiagonalSymmetric{T,L}, B::PairedListDiagonalSymmetric{T,L})
      if A.labels != B.labels || A.nelements != B.nelements
        return($(bin)(full(A), full(B)))
      end
      PairedListDiagonalSymmetric{T,L}($(bin)(A.list, B.list), A.labels, A.nelements)
    end

    function $(bin){T,L}(A::PairedListSymmetric{T,L}, B::PairedListSymmetric{T,L})
      if A.labels != B.labels || A.nelements != B.nelements
        return($(bin)(full(A), full(B)))
      end
      PairedListSymmetric{T,L}($(bin)(A.list, B.list), $(bin)(A.diagonal, B.diagonal), A.labels, A.nelements)
    end

    function $(bin){T,L}(A::PairedListDiagonalSquareTriangular{T,L}, B::PairedListDiagonalSquareTriangular{T,L})
      if A.labels != B.labels || A.nelements != B.nelements || A.upper != B.upper
        return($(bin)(full(A), full(B)))
      end
      PairedListDiagonalSquareTriangular{T,L}($(bin)(A.list, B.list), A.labels, A.nelements, A.upper)
    end

    function $(bin){T,L}(A::PairedListSquareTriangular{T,L}, B::PairedListSquareTriangular{T,L})
      if A.labels != B.labels || A.nelements != B.nelements || A.upper != B.upper
        return($(bin)(full(A), full(B)))
      end
      PairedListSquareTriangular{T,L}($(bin)(A.list, B.list), $(bin)(A.diagonal, B.diagonal), A.labels, A.nelements, A.upper)
    end

  end

end
