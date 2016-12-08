"""
`PairwiseListMatrix{T, diagonal}` is a (squared) symmetric matrix that stores a `list` of
values of type `T` for the pairwise comparison/evaluation of `nelements`.
If `diagonal` is `true` the first element of the list is `1, 1` otherwise is `1, 2`.
If `diagonal` is `false` the diagonal values are stored in a vector on the `diag` field.
"""
type PairwiseListMatrix{T, diagonal} <: AbstractArray{T, 2}
    list::Vector{T}
    diag::Vector{T}
    nelements::Int
end

"Returns `true` if the list has diagonal values."
@inline hasdiagonal{T, diagonal}(::PairwiseListMatrix{T, diagonal}) = diagonal

"Retuns the `list` vector."
@inline getlist(plm::PairwiseListMatrix) = plm.list

"Retuns the `diag` vector (which contains the diagonal values if `diagonal` is `true`)."
@inline getdiag(plm::PairwiseListMatrix) = plm.diag

# listlength
# ==========

"Returns the length of the `list` field"
lengthlist(plm::PairwiseListMatrix) = length(getlist(plm))

function lengthlist{diagonal}(nelements::Int, ::Type{Val{diagonal}})
    diagonal ? div(nelements*(nelements+1),2) : div(nelements*(nelements-1),2)
end

"""
Returns the list length needed for a pairwise measures or comparisons of `nelements`.
If `diagonal` is `true`, diagonal values are included in the list.

```
julia> using PairwiseListMatrices

julia> PLM = PairwiseListMatrix([1, 2, 3, 4, 5, 6], false)
4x4 PairwiseListMatrices.PairwiseListMatrix{Int64,false}:
0  1  2  3
1  0  4  5
2  4  0  6
3  5  6  0

julia> lengthlist(4, false)
6

julia> lengthlist(PLM)
6

```
"""
lengthlist(nelements::Int, diagonal::Bool) = lengthlist(nelements, Val{diagonal})

# Indexes: ij2k
# =============

function ij2k{diagonal}(i::Int, j::Int, nelements::Int, ::Type{Val{diagonal}})
    div(diagonal ?  (nelements*(nelements+1))-((nelements-i)*(nelements-i+1)) :
        (nelements*(nelements-1))-((nelements-i)*(nelements-i-1)), 2) - nelements + j
end

"""
Returns the `k` index of the `list` from the indixes `i` and `j` with `i<j` from a matrix
of `nelements` by `nelements`. `diagonal` should be `true` or `Val{true}` if the diagonal
values are on the `list`. You must not use it with `i>j`.

```
julia> PLM = PairwiseListMatrix([10,20,30,40,50,60], true)
3x3 PairwiseListMatrices.PairwiseListMatrix{Int64,true}:
10  20  30
20  40  50
30  50  60

julia> ij2k(1, 2, 3, true)
2

julia> getlist(PLM)[2]
20

```
"""
ij2k(i::Int, j::Int, nelements::Int, diagonal::Bool) = ij2k(i, j, nelements, Val{diagonal})

# Empties
# -------

"""
An empty `PairwiseListMatrix` can be created for a given `Type` and a number of elements
`nelements`. The `diagonal` (default to `false`) can be declared as `true` to indicate that
the list needs space for the diagonal elements. If `diagonal` is `false` the diagonal
values are stored in a vector on the `diag` field instead of being on the list.
The `diag` vector can be filled with the optional `diagonalvalue` argument (default to `0`).
"""
function (::Type{PairwiseListMatrix}){T}(::Type{T},
                                         nelements::Int,
                                         diagonal::Bool = false,
                                         diagonalvalue::T = zero(T))
    if diagonal
        return PairwiseListMatrix{T, true}(Array(T, lengthlist(nelements, Val{true})),
                                           T[],
                                           nelements)
    else
        return PairwiseListMatrix{T, false}(Array(T, lengthlist(nelements, Val{false})),
                                            fill!(Array(T, nelements), diagonalvalue),
                                            nelements)
    end
end

# Convert
# -------

function Base.convert{T,F,D}(::Type{PairwiseListMatrix{T,D}}, mat::PairwiseListMatrix{F,D})
    PairwiseListMatrix{T,D}(convert(Vector{T}, mat.list),
                            convert(Vector{T}, mat.diag),
                            copy(mat.nelements))
end

function Base.convert{T,F}(::Type{PairwiseListMatrix{T,true}}, mat::PairwiseListMatrix{F,false})
    N = copy(mat.nelements)
    list = Array(T, lengthlist(N, Val{true}))
    k = 0
    @inbounds for i in 1:N
        for j in i:N
             list[k += 1] = mat[i,j]
        end
    end
    PairwiseListMatrix{T, true}(list, T[], N)
end

function Base.convert{T,F}(::Type{PairwiseListMatrix{T, false}}, mat::PairwiseListMatrix{F, true})
    N = copy(mat.nelements)
    diag = Array(T, N)
    matlist = mat.list
    list = Array(T, length(matlist) - N)
    l = 0
    d = 0
    m = 0
    @inbounds for i in 1:N
        for j in i:N
            if i != j
                @inbounds list[l += 1] = matlist[m += 1]
            else
                @inbounds diag[d += 1] = matlist[m += 1]
            end
        end
    end
    PairwiseListMatrix{T, false}(list, diag, N)
end

function Base.convert{T,F,D}(::Type{PairwiseListMatrix{T, D}}, mat::Matrix{F})
    N = size(mat, 1)
    if N != size(mat, 2)
        throw(ErrorException("It should be a squared matrix."))
    end
    plm = PairwiseListMatrix(T, N, D)
    @inbounds for i in 1:N
        for j in i:N
            value = mat[i,j]
            if i != j && value != mat[j,i]
                throw(ErrorException("It should be a symmetric matrix."))
            end
            plm[i,j] = value
        end
    end
    plm
end

# This is faster than list comprehension (2.4 x)
function Base.convert{T,F}(::Type{Array{F,2}}, lm::PairwiseListMatrix{T, true})
    N = lm.nelements
    complete = Array(F, N, N)
    list = lm.list
    k = 0
    @inbounds for col in 1:N
        for row in col:N
            complete[row, col] = complete[col, row] = F(list[k += 1])
        end
    end
    complete
end

function Base.convert{T,F}(::Type{Array{F,2}}, lm::PairwiseListMatrix{T, false})
    N = lm.nelements
    complete = Array(F, N, N)
    list = lm.list
    diag = lm.diag
    k = 0
    l = 0
    @inbounds for col in 1:(N-1)
        for row in (col+1):N
            complete[row, col] = complete[col, row] = F(list[k += 1])
        end
    end
    @inbounds @simd for i in 1:N
        complete[i, i] = F(diag[i])
    end
    complete
end

# full
# ====

"""
Returns a full dense matrix.
This converts a `PairwiseListMatrix{T, D}` into a `Matrix{T}`
"""
Base.full{T,D}(lm::PairwiseListMatrix{T, D}) = convert(Array{T, 2}, lm)

Base.full{T,D}(m::Symmetric{T, PairwiseListMatrix{T,D}}) = full(m.data)

# From a list
# -----------

@inline _nelements(len::Int) = @fastmath div(1+Int(sqrt(1+8*len)),2)
@inline _nelements_with_diagonal(len::Int) = @fastmath div(Int(sqrt(1+8*len)-1),2)

"""
A `PairwiseListMatrix` can be created from a `list`. The `diagonal` (default to `false`)
could be declared as `true` to indicate that the list has the diagonal elements. If
`diagonal` is `false`, the diagonal values are stored in a vector on the `diag` field
instead of being on the list. The `diag` vector can be filled with the optional
`diagonalvalue` argument (default to `0`).
"""
function (::Type{PairwiseListMatrix}){T}(list::AbstractVector{T},
                                         diagonal::Bool = false,
                                         diagonalvalue::T = zero(T))
    if diagonal
        nelements = _nelements_with_diagonal(length(list))
        return PairwiseListMatrix{T, true}(list, T[], nelements)
    else
        nelements = _nelements(length(list))
        return PairwiseListMatrix{T, false}(list,
                                            fill!(Array(T, nelements), diagonalvalue),
                                            nelements)
    end
end

# AbstractArray methods
# =====================

Base.size(m::PairwiseListMatrix) = (m.nelements, m.nelements)
Base.length(m::PairwiseListMatrix) = m.nelements * m.nelements

Base.eltype{T, diagonal}(m::PairwiseListMatrix{T, diagonal}) = T

for F in [:similar, :copy, :zeros, :ones]
    @eval begin
        function Base.$F{T, diagonal}(m::PairwiseListMatrix{T, diagonal})
            PairwiseListMatrix{T, diagonal}($F(m.list),
                                            $F(m.diag),
                                            copy(m.nelements))
        end
    end
end

function Base.similar{T, diagonal, S}(m::PairwiseListMatrix{T, diagonal}, ::Type{S})
    PairwiseListMatrix{S, diagonal}(similar(m.list, S),
                                    similar(m.diag, S),
                                    copy(m.nelements))
end

# TODO:  p = zeros(PairwiseListMatrix{Float64,false}, 3)

# Indexing (getindex)
# ===================

function Base.getindex{T}(lm::PairwiseListMatrix{T, true}, i::Int, j::Int)
    if i <= j
        return(lm.list[ij2k(i, j, lm.nelements, Val{true})])
    else
        return(lm.list[ij2k(j, i, lm.nelements, Val{true})])
    end
end

function Base.getindex{T}(lm::PairwiseListMatrix{T, false}, i::Int, j::Int)
    if i < j
        return(lm.list[ij2k(i, j, lm.nelements, Val{false})])
    elseif i > j
        return(lm.list[ij2k(j, i, lm.nelements, Val{false})])
    else
        return(lm.diag[i])
    end
end

# Base.linearindexing(m::PairwiseListMatrix) = Base.LinearFast()

function Base.getindex{T}(lm::PairwiseListMatrix{T, true}, i::Int)
    n = lm.nelements
    row = Int(ceil(i/n))
    col = i - n*(row-1)
    if row <= col
        @inbounds return lm.list[ij2k(row, col, n, Val{true})]
    else
        @inbounds return lm.list[ij2k(col, row, n, Val{true})]
    end
end

function Base.getindex{T}(lm::PairwiseListMatrix{T, false}, i::Int)
    n = lm.nelements
    row = Int(ceil(i/n))
    col = i - n*(row-1)
    if row < col
        @inbounds return lm.list[ij2k(row, col, n, Val{false})]
    elseif row > col
        @inbounds return lm.list[ij2k(col, row, n, Val{false})]
    else
        @inbounds return lm.diag[row]
    end
end

# Set values (setindex!)
# ======================

function Base.setindex!{T}(lm::PairwiseListMatrix{T, true}, v, i::Int, j::Int)
    if i <= j
        return setindex!(lm.list, v, ij2k(i, j, lm.nelements, Val{true}))
    else
        return setindex!(lm.list, v, ij2k(j, i, lm.nelements, Val{true}))
    end
end

function Base.setindex!{T}(lm::PairwiseListMatrix{T, false}, v, i::Int, j::Int)
    if i < j
        return setindex!(lm.list, v, ij2k(i, j, lm.nelements, Val{false}))
    elseif i > j
        return setindex!(lm.list, v, ij2k(j, i, lm.nelements, Val{false}))
    else
        return setindex!(lm.diag, v, i)
    end
end

function Base.setindex!{T}(lm::PairwiseListMatrix{T, true}, v, i::Int)
    n = lm.nelements
    row = Int(ceil(i/n))
    col = i - n*(row-1)
    if row <= col
        return setindex!(lm.list, v, ij2k(row, col, n, Val{true}))
    else
        return setindex!(lm.list, v, ij2k(col, row, n, Val{true}))
    end
end

function Base.setindex!{T}(lm::PairwiseListMatrix{T, false}, v, i::Int)
    n = lm.nelements
    row = Int(ceil(i/n))
    col = i - n*(row-1)
    if row < col
        return setindex!(lm.list, v, ij2k(row, col, n, Val{false}))
    elseif row > col
        return setindex!(lm.list, v, ij2k(col, row, n, Val{false}))
    else
        return setindex!(lm.diag, v, row)
    end
end

# Transpose
# =========

Base.transpose(lm::PairwiseListMatrix) = lm
Base.transpose!(lm::PairwiseListMatrix) = lm

Base.ctranspose(lm::PairwiseListMatrix) = lm
Base.ctranspose!(lm::PairwiseListMatrix) = lm

# Diagonal
# ========

diagonal{T}(lm::PairwiseListMatrix{T, false}) = lm.diag

diagonal{T}(lm::PairwiseListMatrix{T, true}) = T[ lm[i,i] for i in 1:lm.nelements ]

# Unary operations (faster than default methods)
# ==============================================

for F in (:abs, :-)
    @eval begin
        function Base.$F{T}(lm::PairwiseListMatrix{T, true})
            list = $F(lm.list)
            PairwiseListMatrix{eltype(list), true}(list, lm.diag, lm.nelements)
        end
        function Base.$F{T}(lm::PairwiseListMatrix{T, false})
            list = $F(lm.list)
            diag = $F(lm.diag)
            PairwiseListMatrix{eltype(list), true}(list, diag, lm.nelements)
        end
    end
end

Base.svd(m::PairwiseListMatrix) = svd(full(m))

# Binary operations (faster than default methods)
# ===============================================

for F in ( :-,  :.-, :+, :.+, :.*, :./ )
    @eval begin
        function Base.$F{T,S}(A::PairwiseListMatrix{T, true}, B::PairwiseListMatrix{S, true})
            @assert A.nelements == B.nelements "Matrices must have the same number of elements"
            list = $F(A.list, B.list)
            PairwiseListMatrix{eltype(list), true}(list, copy(A.diag), copy(A.nelements))
        end

        function Base.$F{T,S}(A::PairwiseListMatrix{T, false}, B::PairwiseListMatrix{S, false})
            @assert A.nelements == B.nelements "Matrices must have the same number of elements"
            list = $F(A.list, B.list)
            diag = $F(A.diag, B.diag)
            PairwiseListMatrix{eltype(list), false}(list, diag, copy(A.nelements))
        end
    end
end

for F in ( :-,  :.-, :+, :.+, :.*, :*, :/, :./ )
    @eval begin
        function Base.$F{T, S <: Number}(A::PairwiseListMatrix{T, true}, B::S)
            list = $F(A.list, B)
            PairwiseListMatrix{eltype(list), true}(list, copy(A.diag), copy(A.nelements))
        end

        function Base.$F{T, S <: Number}(A::PairwiseListMatrix{T, false}, B::S)
            list = $F(A.list, B)
            diag = $F(A.diag, B)
            PairwiseListMatrix{eltype(list), false}(list, diag, copy(A.nelements))
        end

        Base.$F{T,D,S <: Number}(B::S, A::PairwiseListMatrix{T,D}) = $F(A, B)
    end
end

# Faster mean and sum
# ===================

Base.sum{T}(m::PairwiseListMatrix{T, false}) = T(2)*sum(m.list) + sum(m.diag)
Base.sum{T}(m::PairwiseListMatrix{T, true}) =  T(2)*sum(m.list) - sum(diagonal(m))

function _sum_kernel!(sum_i, list, N)
    k = 0
    l = 0
    for i in 1:N
        l += 1
        @inbounds for j in i:N
            sum_i[i] += list[k += 1]
        end
        @inbounds for j in (i+1):N
            sum_i[j] += list[l += 1]
        end
    end
    sum_i
end

function _sum_kernel!(sum_i, diag, list, N)
    k = 0
    l = 0
    for i in 1:(N-1)
        @inbounds for j in (i+1):N
            sum_i[i] += list[k += 1]
        end
        @inbounds for j in (i+1):N
            sum_i[j] += list[l += 1]
        end
    end
    @inbounds @simd for i in 1:N
        sum_i[i] += diag[i]
    end
    sum_i
end

function _test_and_set_sum{T}(::Type{T}, region, N)
    if region == 1
        sum_i = zeros(T, 1, N)
    elseif region == 2
        sum_i = zeros(T, N, 1)
    else
        throw(ErrorException("region should be 1 or 2"))
    end
    sum_i
end

function Base.sum{T}(lm::PairwiseListMatrix{T, true}, region::Int)
    N = lm.nelements
    sum_i = _test_and_set_sum(T, region, N)
    _sum_kernel!(sum_i, lm.list, N)
end

function Base.sum{T}(lm::PairwiseListMatrix{T, false}, region::Int)
    N = lm.nelements
    sum_i = _test_and_set_sum(T, region, N)
    _sum_kernel!(sum_i, lm.diag, lm.list, N)
end

Base.mean(m::PairwiseListMatrix) = sum(m)/length(m)
Base.mean(m::PairwiseListMatrix, region::Int) = sum(m, region) ./ m.nelements

# Sum/Mean without diagonal: sum/mean_nodiag
# ------------------------------------------

"Sum the values outside the diagonal"
sum_nodiag{T}(m::PairwiseListMatrix{T, false}) = T(2) * sum(m.list)
sum_nodiag{T}(m::PairwiseListMatrix{T, true}) = T(2) * sum(m.list) - sum(diagonal(m))

function _sum_nodiag_kernel!{T}(sum_i, lm::PairwiseListMatrix{T, true}, N)
    list = lm.list
    k = 0
    @inbounds for i in 1:N
        for j in i:N
            k += 1
            if i != j
                value = list[k]
                sum_i[i] += value
                sum_i[j] += value
            end
        end
    end
    sum_i
end

function _sum_nodiag_kernel!{T}(sum_i, lm::PairwiseListMatrix{T, false}, N)
    list = lm.list
    k = 0
    l = 0
    for i in 1:(N-1)
        @inbounds for j in (i+1):N
            sum_i[i] += list[k += 1]
        end
        @inbounds for j in (i+1):N
            sum_i[j] += list[l += 1]
        end
    end
    sum_i
end

function sum_nodiag{T, diagonal}(lm::PairwiseListMatrix{T, diagonal}, region::Int)
    N = lm.nelements
    sum_i = _test_and_set_sum(T, region, N)
    _sum_nodiag_kernel!(sum_i, lm, N)
end

"Mean of the values outside the diagonal"
mean_nodiag(m::PairwiseListMatrix) = sum_nodiag(m) / (length(m) - m.nelements)
mean_nodiag(m::PairwiseListMatrix, region::Int) = sum_nodiag(m, region) ./ (m.nelements-1)

# Operations on Vector{PairwiseListMatrix}
# ========================================

# sum
# ---

@inline _has_diagonal{T, diagonal}(x::PairwiseListMatrix{T, diagonal}) = diagonal

function Base.sum{T, diagonal}(list::Vector{PairwiseListMatrix{T, diagonal}})
    samples = length(list)
    if samples == 0
        throw(ErrorException("Empty list"))
    end
    if samples == 1
        return list[1]
    else
        start = copy(list[1])
        i = 2
        start_list = start.list
        N = length(start_list)
        i = 2
        while i <= samples
            @inbounds ylist = list[i].list
            if length(ylist) != N
                throw(ErrorException("Different number of elements"))
            end
            @inbounds @simd for k in 1:N
                start_list[k] += ylist[k]
            end
            i += 1
        end
        if !diagonal
            i = 2
            start_diag = start.diag
            while i <= samples
                @inbounds ydiag = list[i].diag
                @inbounds @simd for k in 1:length(start_diag)
                    start_diag[k] += ydiag[k]
                end
                i += 1
            end
        end
        return start
    end
end

# std
# ---

function Base.varm{T, D}(list::Vector{PairwiseListMatrix{T, D}}, mean::PairwiseListMatrix{T, D})
    samples = length(list)
    if samples < 2
        throw(ErrorException("You need at least 2 samples."))
    end
    out = zeros(list[1])
    out_list = out.list
    mean_list = mean.list
    N = length(mean_list)
    @inbounds for sample in list
        sample_list = sample.list
        if length(sample_list) != N
            throw(ErrorException("Different number of elements"))
        end
        @inbounds @simd for k in 1:N
            out_list[k] += abs2( sample_list[k] - mean_list[k])
        end
    end
    if !D
        out_diag = out.diag
        mean_diag = mean.diag
        @inbounds for sample in list
            sample_diag = sample.diag
            for k in 1:length(out_diag)
                @inbounds out_diag[k] += abs2( sample_diag[k] - mean_diag[k])
            end
        end
    end
    out ./ samples
end

function Base.var{T, D}(list::Vector{PairwiseListMatrix{T, D}}; mean=nothing)
    if mean === nothing
        varm(list, Base.mean(list))
    else
        varm(list, mean)
    end
end

function Base.std{T, D}(list::Vector{PairwiseListMatrix{T, D}}; mean=nothing)
    sqrt(var(list, mean=mean))
end

# zscore
# ======

function zscore!{E <: AbstractFloat, D}(list::Vector{PairwiseListMatrix{E,D}},
                                        mat::PairwiseListMatrix{E,D})
    list_mean = mean(list)
    if size(mat) != size(list_mean)
        throw(ErrorException("PairwiseListMatrices must have the same size"))
    end
    list_std  = std(list, mean=list_mean)

    mat_list = mat.list
    std_list = list_std.list
    mean_list = list_mean.list

    ElementType = eltype(mat_list)

    @inbounds for i in 1:length(mat_list)
        s = std_list[i]
        m = mean_list[i]
        value = mat_list[i]
        if (m - value + one(ElementType)) ≉ one(ElementType)
            if (s + one(ElementType)) ≉ one(ElementType)
                mat_list[i] = (value - m)/s
            else
                mat_list[i] = NaN
            end
        else
            mat_list[i] = zero(ElementType)
        end
    end

    if D
        mat_diag = mat.diag
        std_diag = list_std.diag
        mean_diag = list_mean.diag
        @inbounds for i in 1:length(mat_diag)
            s = std_diag[i]
            m = mean_diag[i]
            value = mat_diag[i]
            if (m - value + one(ElementType)) ≉ one(ElementType)
                if (s + one(ElementType)) ≉ one(ElementType)
                    mat_diag[i] = (value - m)/s
                else
                    mat_diag[i] = NaN
                end
            else
                mat_diag[i] = zero(ElementType)
            end
        end
    end

    mat
end

zscore(list, mat) = zscore!(list, copy(mat))

# Tables
# ======

"""
Creates a `Matrix{Any}` is useful for `writedlm` and/or `writecsv`.
Labels are stored in the columns 1 and 2, and the values in the column 3.
Diagonal values are included by default.

```
julia> list
3x3 PairwiseListMatrices.PairwiseListMatrix{Int64,false}:
0  1  2
1  0  3
2  3  0

julia> to_table(list)
6x3 Array{Any,2}:
"A"  "A"  0
"A"  "B"  1
"A"  "C"  2
"B"  "B"  0
"B"  "C"  3
"C"  "C"  0

julia> to_table(list, false)
3x3 Array{Any,2}:
"A"  "B"  1
"A"  "C"  2
"B"  "C"  3

```
"""
function to_table(plm::PairwiseListMatrix, diagonal::Bool=true)
    N = plm.nelements
    labels = 1:N
    table = Array(Any, diagonal ? div(N*(N+1),2) : div(N*(N-1),2), 3)
    t = 0
    @iterateupper plm diagonal begin
        t += 1
        table[t, 1] = labels[i]
        table[t, 2] = labels[j]
        table[t, 3] = list[k]
    end
    table
end

"""
Creation of a `PairwiseListMatrix` from a `Matrix`.
By default the columns with the labels for i (slow) and j (fast) are 1 and 2.
Values are taken from the column 3 by default.

```
julia> data = readcsv("example.csv")
3x3 Array{Any,2}:
"A"  "B"  10
"A"  "C"  20
"B"  "C"  30

julia> from_table(data, Int, false)
3x3 PairwiseListMatrices.PairwiseListMatrix{Int64,false}:
0  10  20
10   0  30
20  30   0

```
"""
function from_table{T}(table::Matrix,
                       value::Type{T},
                       diagonal::Bool,
                       labelcols::Vector{Int}=[1,2],
                       valuecol::Int=3)
    PairwiseListMatrix(convert(Vector{T}, table[:,valuecol]),
                       # unique(table[:,labelcols]),
                       diagonal)
end

# write...
# ========

"""
This function takes the filename as first argument and a `PairwiseListMatrix` as second argument.
If the third positional argument is `true` (default) the diagonal is included in the output.
The keyword argument `delim` (by default is `'\t’`) allows to modified the character used as delimiter.

julia> PLM  = PairwiseListMatrix(collect(1:6), true)
3x3 PairwiseListMatrices.PairwiseListMatrix{Int64,true}:
1  2  3
2  4  5
3  5  6

julia> writedlm("example.txt", PLM, false, delim=' ')

shell> cat example.txt
1 2 2
1 3 3
2 3 5

"""
function Base.writedlm(filename::String,
                  plm::PairwiseListMatrix,
                  diagonal::Bool=true;
                  delim::Char='\t')
    open(filename, "w") do fh
        # if length(plm.labels) != 0
        #    labels = plm.labels
        #    @iterateupper plm diagonal println(fh, labels[i], delim, labels[j], delim, list[k])
        # else
        @iterateupper plm diagonal println(fh, i, delim, j, delim, list[k])
        # end
    end
end

"""
This function takes the filename as first argument and a `PairwiseListMatrix` as second argument.
If the third positional argument is `true` (default) the diagonal is included in the output.

julia> PLM  = PairwiseListMatrix(collect(1:6), true)
3x3 PairwiseListMatrices.PairwiseListMatrix{Int64,true}:
1  2  3
2  4  5
3  5  6

julia> writecsv("example.csv", PLM)

shell> cat example.csv
1,1,1
1,2,2
1,3,3
2,2,4
2,3,5
3,3,6

julia> writecsv("example.csv", PLM, false)

shell> cat example.csv
1,2,2
1,3,3
2,3,5

"""
function Base.writecsv(filename::String, plm::PairwiseListMatrix, diagonal::Bool=true)
    writedlm(filename, plm, diagonal, delim=',')
end

# triu! & triu
# ============

function Base.triu!(mat::PairwiseListMatrix)
    throw(ErrorException("PairwiseListMatrix must be Symmetric, use triu instead of triu!"))
end
Base.triu!(mat::PairwiseListMatrix, k::Int) = triu!(mat)
Base.triu(mat::PairwiseListMatrix) = triu(full(mat))
Base.triu(mat::PairwiseListMatrix, k::Int) = triu(full(mat), k::Int)

# # join
# # ====
#
# """
# This function join two PairwiseListMatrices by their labels, returning two PairwiseListMatrices with same size and labels.
# There are 4 `kind` of joins:
# - `:inner` : Intersect. The output matrices only include the labels that are in both PairwiseListMatrices
# - `:outer` : Union. Include the labels of the two PairwiseListMatrices.
# - `:left` : Only use labels form the first argument.
# - `:right` : Only use labels form the second argument.
# `NaN`s are filled in where needed to complete joins. The default value for missing values can be changed passing a tuple to `missing`.
#
# ```
# julia> l = PairwiseListMatrix([1.,2.,3.], ['a', 'b', 'c'], false) # a b c
# 3x3 PairwiseListMatrices.PairwiseListMatrix{Float64,false}:
# 0.0  1.0  2.0
# 1.0  0.0  3.0
# 2.0  3.0  0.0
#
# julia> r = PairwiseListMatrix([1.,2.,3.], ['b', 'c', 'd'], false) #   b c d
# 3x3 PairwiseListMatrices.PairwiseListMatrix{Float64,false}:
# 0.0  1.0  2.0
# 1.0  0.0  3.0
# 2.0  3.0  0.0
#
# julia> join(l, r, kind=:inner)                                    #   b c
# (
# 2x2 PairwiseListMatrices.PairwiseListMatrix{Float64,false}:
# 0.0  3.0
# 3.0  0.0,
#
# 2x2 PairwiseListMatrices.PairwiseListMatrix{Float64,false}:
# 0.0  1.0
# 1.0  0.0)
#
# julia> join(l, r, kind=:outer)                                    # a b c d
# (
#     4x4 PairwiseListMatrices.PairwiseListMatrix{Float64,false}:
#     0.0    1.0    2.0  NaN
#     1.0    0.0    3.0  NaN
#     2.0    3.0    0.0  NaN
#     NaN    NaN    NaN    NaN,
#
#     4x4 PairwiseListMatrices.PairwiseListMatrix{Float64,false}:
#     NaN  NaN    NaN    NaN
#     NaN    0.0    1.0    2.0
#     NaN    1.0    0.0    3.0
#     NaN    2.0    3.0    0.0)
#
#     ```
#     """
#     function join{L <: AbstractFloat, R <: AbstractFloat, DL, DR}(left::PairwiseListMatrix{L, DL}, right::PairwiseListMatrix{R, DR}; kind::Symbol = :inner ,
#         missing::Tuple{L,R} = (L(NaN), R(NaN)) )
#
#         labels_left  = labels(left)
#         labels_right = labels(right)
#
#         if labels_left == labels_right
#             return(left, right)
#         end
#
#         if     kind == :inner
#             out_labels = intersect(labels_left, labels_right)
#             N = length(out_labels)
#             out_L = PairwiseListMatrix(L, N, out_labels, DL)
#             out_R = PairwiseListMatrix(R, N, out_labels, DR)
#             for i in 1:N
#                 li = out_labels[i]
#                 for j in i:N
#                     lj = out_labels[j]
#                     out_L[i,j] = getlabel(left,  li, lj)
#                     out_R[i,j] = getlabel(right, li, lj)
#                 end
#             end
#             return(out_L, out_R)
#         elseif kind == :left
#             out_labels = labels_left
#             N = length(out_labels)
#             out_R = PairwiseListMatrix(R, N, out_labels, DR)
#             for i in 1:N
#                 li = out_labels[i]
#                 flag_i_r = li in labels_right
#                 for j in i:N
#                     lj = out_labels[j]
#                     out_R[i,j] = (flag_i_r && (lj in labels_right)) ? getlabel(right, li, lj) : missing[2]
#                 end
#             end
#             return(left, out_R)
#         elseif kind == :right
#             out_labels = labels_right
#             N = length(out_labels)
#             out_L = PairwiseListMatrix(L, N, out_labels, DL)
#             for i in 1:N
#                 li = out_labels[i]
#                 flag_i_l = li in labels_left
#                 for j in i:N
#                     lj = out_labels[j]
#                     out_L[i,j] = (flag_i_l  && (lj in labels_left)) ? getlabel(left,  li, lj) : missing[1]
#                 end
#             end
#             return(out_L, right)
#         elseif kind == :outer
#             out_labels = union(labels_left, labels_right)
#             N = length(out_labels)
#             out_L = PairwiseListMatrix(L, N, out_labels, DL)
#             out_R = PairwiseListMatrix(R, N, out_labels, DR)
#             for i in 1:N
#                 li = out_labels[i]
#                 flag_i_l = li in labels_left
#                 flag_i_r = li in labels_right
#                 for j in i:N
#                     lj = out_labels[j]
#                     out_L[i,j] = (flag_i_l && (lj in labels_left))  ? getlabel(left,  li, lj) : missing[1]
#                     out_R[i,j] = (flag_i_r && (lj in labels_right)) ? getlabel(right, li, lj) : missing[2]
#                 end
#             end
#             return(out_L, out_R)
#         else
#             throw(ArgumentError("Unknown kind of join requested: use :inner, :left, :right or :outer"))
#         end
#
#     end
