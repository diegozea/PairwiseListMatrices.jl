"""
`PairwiseListMatrix{T, L, diagonal}` is a (squared) symmetric matrix that stores a `list` of values of type `T` for the pairwise comparison/evaluation of `nelements`.
If `diagonal` is `true` the first element of the list is `1, 1` otherwise is `1, 2`.
If `diagonal` is `false`, the diagonal values are stored in a vector on the `diag` field.
Labels can be stored on the field `labels` as an `IndexedArray`.
"""
type PairwiseListMatrix{T, L, diagonal} <: AbstractArray{T, 2}
  list::Vector{T}
  diag::Vector{T}
  labels::IndexedArray{L}
  nelements::Int
end

# Creation
# ========

function _test_nelements(vector, nelements::Int, name::ASCIIString)
  if length(vector) != 0 && length(vector) != nelements
    throw(DimensionMismatch(string(name, " must have ", nelements, " names!")))
  end
end

@inline _list_length(nelements) = @fastmath div(nelements*(nelements-1),2)
@inline _list_with_diagonal_length(nelements) = @fastmath div(nelements*(nelements+1),2)

# Empties
# -------

"""
An empty `PairwiseListMatrix` can be created for a given `Type` and a number of elements `elements`.
Optionally, you can use a vector or IndexedArray for declaring names/labels to each element.
The `diagonal` (default to `false`) could be declared as `true` in order to indicate that the list needs space for the diagonal elements.
If `diagonal` is `false`, the diagonal values are stored in a vector on the `diag` field instead of being on the list.
This vector can be filled to a value with the optional `diagonalvalue` arguments (default to `0`).
"""
function PairwiseListMatrix{T, L}(::Type{T}, nelements::Int, labels::IndexedArray{L},
                                  diagonal::Bool=false, diagonalvalue::T=zero(T))
  _test_nelements(labels, nelements, "labels")
  if diagonal
    return PairwiseListMatrix{T, L, diagonal}(Array(T, _list_with_diagonal_length(nelements)), T[], labels, nelements)
  else
    return PairwiseListMatrix{T, L, diagonal}(Array(T, _list_length(nelements)), fill!(Array(T, nelements), diagonalvalue), labels, nelements)
  end
end

PairwiseListMatrix{T, L}(::Type{T}, nelements::Int, labels::AbstractVector{L}, diagonal::Bool=false,
                         diagonalvalue::T=zero(T)) = PairwiseListMatrix(T, nelements, IndexedArray(labels), diagonal, diagonalvalue)

PairwiseListMatrix{T}(::Type{T}, nelements::Int, diagonal::Bool=false,
                      diagonalvalue::T=zero(T)) = PairwiseListMatrix(T, nelements, IndexedArray{Any}(), diagonal, diagonalvalue)

# From a list
# -----------

@inline _nelements(len::Int) = @fastmath div(1+Int(sqrt(1+8*len)),2)
@inline _nelements_with_diagonal(len::Int) = @fastmath div(Int(sqrt(1+8*len)-1),2)

"""
A `PairwiseListMatrix` can be from a `list`. Optionally, you can use a vector or IndexedArray for declaring `labels` to each element.
The `diagonal` (default to `false`) could be declared as `true` in order to indicate that the list has the diagonal elements.
If `diagonal` is `false`, the diagonal values are stored in a vector on the `diag` field instead of being on the list.
This vector can be filled to a value with the optional `diagonalvalue` arguments (default to `0`).
"""
function PairwiseListMatrix{T, L}(list::AbstractVector{T}, labels::IndexedArray{L},
                                  diagonal::Bool=false, diagonalvalue::T=zero(T))
  if diagonal
    nelements = _nelements_with_diagonal(length(list))
    _test_nelements(labels, nelements, "labels")
    return PairwiseListMatrix{T, L, diagonal}(list, T[], labels, nelements)
  else
    nelements = _nelements(length(list))
    _test_nelements(labels, nelements, "labels")
    return PairwiseListMatrix{T, L, diagonal}(list, fill!(Array(T, nelements), diagonalvalue), labels, nelements)
  end
end

PairwiseListMatrix{T, L}(list::AbstractVector{T}, labels::AbstractVector{L}, diagonal::Bool=false,
                         diagonalvalue::T=zero(T)) = PairwiseListMatrix(list, IndexedArray(labels), diagonal, diagonalvalue)

PairwiseListMatrix{T}(list::AbstractVector{T}, diagonal::Bool=false,
                      diagonalvalue::T=zero(T)) = PairwiseListMatrix(list, IndexedArray{Any}(), diagonal, diagonalvalue)


# AbstractArray methods
# =====================

size(m::PairwiseListMatrix) = (m.nelements, m.nelements)
length(m::PairwiseListMatrix) = m.nelements * m.nelements

eltype{T, L, diagonal}(m::PairwiseListMatrix{T, L, diagonal}) = T

similar{T, L, diagonal}(m::PairwiseListMatrix{T, L, diagonal}) = PairwiseListMatrix{T, L, diagonal}(similar(m.list), copy(m.diag),
                                                                                                    copy(m.labels), copy(m.nelements))
similar{T, L, diagonal, S}(m::PairwiseListMatrix{T, L, diagonal}, ::Type{S}) = PairwiseListMatrix{S, L, diagonal}(similar(m.list, S), convert(Vector{S}, m.diag),
                                                                                                                  copy(m.labels), S(copy(m.nelements)))

copy{T, L, diagonal}(m::PairwiseListMatrix{T, L, diagonal}) = PairwiseListMatrix{T, L, diagonal}(copy(m.list), copy(m.diag), copy(m.labels), copy(m.nelements))

# Indexing (getindex)
# ===================

@inline _listindex(i, j, n) = @fastmath div((n*(n-1))-((n-i)*(n-i-1)),2) - n + j
@inline _listindex_with_diagonal(i, j, n) = @fastmath div((n*(n+1))-((n-i)*(n-i+1)),2) - n + j

function getindex{T, L}(lm::PairwiseListMatrix{T, L, true}, i::Int, j::Int)
  if i <= j
    return(lm.list[_listindex_with_diagonal(i, j, lm.nelements)])
  else
    return(lm.list[_listindex_with_diagonal(j, i, lm.nelements)])
  end
end

function getindex{T, L}(lm::PairwiseListMatrix{T, L, false}, i::Int, j::Int)
  if i < j
    return(lm.list[_listindex(i, j, lm.nelements)])
  elseif i > j
    return(lm.list[_listindex(j, i, lm.nelements)])
  else
    return(lm.diag[i])
  end
end

# Labels
# ======

"Get the labels/names of the row/columns of the matrix"
labels(lm::PairwiseListMatrix) = lm.labels

"You can use labels for add labels/names to the row/columns of the matrix"
function labels!(lm::PairwiseListMatrix, labels::IndexedArray)
  _test_nelements(labels, lm.nelements, "labels")
  lm.labels = labels
end

labels!(lm::PairwiseListMatrix, labels::AbstractVector) = labels!(lm, IndexedArray(labels))

# Indexing using labels (getlabel)
# ================================

"Like `getindex`, but using the labels/names instead of `Int` numbers."
function getlabel(lm::PairwiseListMatrix, i, j)
  if isempty(lm.labels)
    throw(ErrorException("There are not labels in the matrix. You can use labels!(...) for add them."))
  else
    return getindex(lm, findfirst(lm.labels, i), findfirst(lm.labels, j))
  end
end

# Set values (setindex!)
# ======================

function setindex!{T, L}(lm::PairwiseListMatrix{T, L, true}, v, i::Int, j::Int)
  if i <= j
    return setindex!(lm.list, v, _listindex_with_diagonal(i, j, lm.nelements))
  else
    return setindex!(lm.list, v, _listindex_with_diagonal(j, i, lm.nelements))
  end
end

function setindex!{T, L}(lm::PairwiseListMatrix{T, L, false}, v, i::Int, j::Int)
  if i < j
    return setindex!(lm.list, v, _listindex(i, j, lm.nelements))
  elseif i > j
    return setindex!(lm.list, v, _listindex(j, i, lm.nelements))
  else
    return setindex!(lm.diag, v, i)
  end
end

# Set values using labels (setlabel!)
# ===================================

"Like `setindex!`, but using the labels/names instead of `Int` numbers."
function setlabel!(lm::PairwiseListMatrix, v, i, j)
  if isempty(lm.labels)
    throw(ErrorException("There are not labels in the matrix. You can use labels!(...) for add them."))
  else
    return setindex!(lm, v, findfirst(lm.labels, i), findfirst(lm.labels, j))
  end
end

# Transpose
# =========

transpose(lm::PairwiseListMatrix) = lm
transpose!(lm::PairwiseListMatrix) = lm

ctranspose(lm::PairwiseListMatrix) = lm
ctranspose!(lm::PairwiseListMatrix) = lm

# diag and full
# =============

diag{T, L}(lm::PairwiseListMatrix{T, L, false}) = lm.diag

diag{T, L}(lm::PairwiseListMatrix{T, L, true}) = T[ lm[i,i] for i in 1:lm.nelements ]

# This is faster than list comprehension (1.21 x)
"Returns a full dense matrix"
function full{T, L}(lm::PairwiseListMatrix{T, L, true})
  complete = Array(T, lm.nelements, lm.nelements)
  k = 0
  @inbounds for col in 1:lm.nelements
    for row in col:lm.nelements
      k += 1
      complete[row, col] = lm.list[k]
    end
  end
  @inbounds for col in 2:lm.nelements
    for row in 1:(col-1)
      complete[row, col] = lm[row, col]
    end
  end
  complete
end

function full{T, L}(lm::PairwiseListMatrix{T, L, false})
  complete = Array(T, lm.nelements, lm.nelements)
  k = 0
  @inbounds for col in 1:(lm.nelements-1)
    for row in (col+1):lm.nelements
      k += 1
      complete[row, col] = lm.list[k]
    end
  end
  @inbounds for col in 1:lm.nelements
    for row in 1:col
      complete[row, col] = lm[row, col]
    end
  end
  complete
end

full{T, S <: PairwiseListMatrix}(m::Symmetric{T, S}) = full(m.data)

# Unary operations
# ================

for una in (:abs, :-)
  @eval begin
    $(una){T, L}(lm::PairwiseListMatrix{T, L, true}) = PairwiseListMatrix{T, L, true }($(una)(lm.list), lm.diag, lm.labels, lm.nelements)
    $(una){T, L}(lm::PairwiseListMatrix{T, L, false})= PairwiseListMatrix{T, L, false}($(una)(lm.list), $(una)(lm.diag), lm.labels, lm.nelements)
  end
end

svd(m::PairwiseListMatrix) = svd(full(m))

# Binary operations
# =================

for bin in ( :-, :+, :.*, :./, :.+, :.- )

  @eval begin

    function $(bin){T, L}(A::PairwiseListMatrix{T, L, true}, B::PairwiseListMatrix{T, L, true})
      if A.labels != B.labels || A.nelements != B.nelements
        return($(bin)(full(A), full(B)))
      end
      PairwiseListMatrix{T, L, true}($(bin)(A.list, B.list), A.diag, A.labels, A.nelements)
    end

    function $(bin){T, L}(A::PairwiseListMatrix{T, L, false}, B::PairwiseListMatrix{T, L, false})
      if A.labels != B.labels || A.nelements != B.nelements
        return($(bin)(full(A), full(B)))
      end
      PairwiseListMatrix{T, L, false}($(bin)(A.list, B.list), $(bin)(A.diag, B.diag), A.labels, A.nelements)
    end

    $(bin)(A::PairwiseListMatrix, B::PairwiseListMatrix) = $(bin)(full(A), full(B))

  end

end

for bin in (:*, :/)

  @eval $(bin)(A::PairwiseListMatrix, B::PairwiseListMatrix) = $(bin)(full(A), full(B))

end

# Fast operations
# ===============

mean{T, L}(m::PairwiseListMatrix{T, L, false}) = (2*sum(m.list) + sum(m.diag))/length(m)
