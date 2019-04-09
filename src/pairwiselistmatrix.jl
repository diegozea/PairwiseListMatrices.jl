"""
`PairwiseListMatrix{T, diagonal, VT}` is a (squared) symmetric matrix that stores a `list`
of type `VT` with values of type `T` for the pairwise comparison/evaluation of `nelements`.
If `diagonal` is `true` the first element of the list is `1, 1` otherwise is `1, 2`.
If `diagonal` is `false` the diagonal values are stored in a vector on the `diag` field.
"""
mutable struct PairwiseListMatrix{T,diagonal,VT} <: AbstractArray{T, 2}
    list::VT
    diag::VT
    nelements::Int

    function PairwiseListMatrix{T,diagonal,VT}(list::VT,
                                               diag::VT,
                                               N::Int) where {T,diagonal,VT}
        @assert eltype(VT) === T "Field vectors must have the same element type than the PairwiseListMatrix"
        if !diagonal
            n_diag = length(diag)
            n_list = _nelements(length(list))
            @assert n_diag == N "There are $n_diag elements in the diag field instead of $N"
        else
            n_list = _nelements_with_diagonal(length(list))
            @assert n_list == N "There are $n_list elements in the list field instead of $N"
        end
        new(list, diag, N)
    end
end

"Returns `true` if the list has diagonal values."
@inline hasdiagonal(::PairwiseListMatrix{T,diagonal,VT}) where {T,diagonal,VT} = diagonal
@inline hasdiagonal(nplm::NamedArray{T,2,PairwiseListMatrix{T,D,TV},DN}) where {T,D,TV,DN} = D

"Retuns the `list` vector."
@inline getlist(plm::PairwiseListMatrix) = plm.list

"Retuns the `diag` vector (which contains the diagonal values if `diagonal` is `false`)."
@inline getdiag(plm::PairwiseListMatrix) = plm.diag

# listlength
# ==========

"Returns the length of the `list` field"
lengthlist(plm::PairwiseListMatrix) = length(getlist(plm))

function lengthlist(nelements::Int, ::Type{Val{diagonal}}) where diagonal
    diagonal ? div(nelements*(nelements+1),2) : div(nelements*(nelements-1),2)
end

"""
Returns the list length needed for a pairwise measures or comparisons of `nelements`.
If `diagonal` is `true`, diagonal values are included in the list.

```jldoctest
julia> using PairwiseListMatrices

julia> plm = PairwiseListMatrix([1, 2, 3, 4, 5, 6], false)
4×4 PairwiseListMatrix{Int64,false,Array{Int64,1}}:
 0  1  2  3
 1  0  4  5
 2  4  0  6
 3  5  6  0

julia> lengthlist(4, false)
6

julia> lengthlist(plm)
6

```
"""
lengthlist(nelements::Int, diagonal::Bool) = lengthlist(nelements, Val{diagonal})

# Indexes: ij2k
# =============

@inline @fastmath function ij2k(i::Int, j::Int, nelements::Int,
                                ::Type{Val{diagonal}}) where diagonal
    div(diagonal ?  (nelements*(nelements+1))-((nelements-i)*(nelements-i+1)) :
    (nelements*(nelements-1))-((nelements-i)*(nelements-i-1)), 2) - nelements + j
end

"""
Returns the `k` index of the `list` from the indixes `i` and `j` with `i<j` from a matrix
of `nelements` by `nelements`. `diagonal` should be `true` or `Val{true}` if the diagonal
values are on the `list`. You must not use it with `i>j`.

```jldoctest
julia> using PairwiseListMatrices

julia> plm = PairwiseListMatrix([10,20,30,40,50,60], true)
3×3 PairwiseListMatrix{Int64,true,Array{Int64,1}}:
 10  20  30
 20  40  50
 30  50  60

julia> ij2k(1, 2, 3, true)
2

julia> getlist(plm)[2]
20

```
"""
@inline function ij2k(i::Int, j::Int, nelements::Int, diagonal::Bool)
    ij2k(i, j, nelements, Val{diagonal})
end

# Empties
# -------

function _diagonal_value(diagonal::Bool, ::Type{T}) where T
    if !diagonal
        if hasmethod(zero, (T,))
            return zero(T)
        elseif T <: Union{Missing, Int}
            return missing
        elseif T === Any
            return nothing
        else
            throw(ErrorException(
            "Please use the last argument to fill the diagonal. It should be of type $T."))
        end
    end
    Array{T}(undef, 1)[1]
end


"""
An empty `PairwiseListMatrix` can be created for a given `Type` and a number of elements
`nelements`. The `diagonal` (default to `false`) can be declared as `true` to indicate that
the list needs space for the diagonal elements. If `diagonal` is `false` the diagonal
values are stored in a vector on the `diag` field instead of being on the list.
The `diag` vector can be filled with the optional `diagonalvalue` argument (default to `0`).

```julia
PairwiseListMatrix(Int, 3)
PairwiseListMatrix(Int, 3, true)
```
"""
function PairwiseListMatrix(::Type{T},
                           nelements::Int,
                           diagonal::Bool = false,
                           diagonalvalue::T = _diagonal_value(diagonal, T)) where T
    if diagonal
        return PairwiseListMatrix{T, true, Vector{T}}(
            Array{T}(undef, lengthlist(nelements, Val{true})),
            T[],
            nelements)
    else
        return PairwiseListMatrix{T, false, Vector{T}}(
            Array{T}(undef, lengthlist(nelements, Val{false})),
            fill!(Array{T}(undef, nelements), diagonalvalue),
            nelements)
    end
end

"""
```julia
PairwiseListMatrix{Int,false,Vector{Int}}(3)
```
"""
function PairwiseListMatrix{T,diagonal,VT}(nelements::Int) where {T,diagonal,VT}
    if diagonal
        return PairwiseListMatrix{T,true,VT}(
            convert(VT, Array{T}(undef, lengthlist(nelements, Val{true}))),
            convert(VT, T[]),
            nelements)
    else
        return PairwiseListMatrix{T, false, VT}(
            Array{T}(undef, lengthlist(nelements, Val{false})),
            convert(VT, Array{T}(undef, nelements)),
            nelements)
    end
end

# Convert
# -------

function Base.convert(::Type{PairwiseListMatrix{T,D,VT}},
                     mat::PairwiseListMatrix{F,D,VF}) where {T,F,D,VT,VF}
    PairwiseListMatrix{T,D,VT}(convert(VT, mat.list),
        convert(VT, mat.diag),
        copy(mat.nelements))
end

function Base.convert(::Type{PairwiseListMatrix{T,true,VT}},
                     mat::PairwiseListMatrix{F,false,VF}) where {T,F,VT,VF}
    N = copy(mat.nelements)
    list = convert(VT, Array{T}(undef, lengthlist(N, Val{true})))
    k = 0
    @inbounds for i in 1:N
        for j in i:N
            list[k += 1] = mat[i,j]
        end
    end
    PairwiseListMatrix{T, true, VT}(list, convert(VT,T[]), N)
end

function Base.convert(::Type{PairwiseListMatrix{T, false, VT}},
                     mat::PairwiseListMatrix{F, true, VF}) where {T,F,VT,VF}
    N = copy(mat.nelements)
    diag = convert(VT, Array{T}(undef, N))
    matlist = mat.list
    list = convert(VT, Array{T}(undef, length(matlist) - N))
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
    PairwiseListMatrix{T, false, VT}(list, diag, N)
end

function Base.convert(::Type{PairwiseListMatrix{T,D,VT}}, mat::Matrix{F}) where {T,F,D,VT}
    N = size(mat, 1)
    if N != size(mat, 2)
        throw(ErrorException("It should be a squared matrix."))
    end
    plm = PairwiseListMatrix{T,D,VT}(N)
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
function Base.convert(::Type{Array{F,2}}, lm::PairwiseListMatrix{T, true, VT}) where {T,F,VT}
    N = lm.nelements
    complete = Array{F}(undef, N, N)
    list = lm.list
    k = 0
    @inbounds for col in 1:N
        for row in col:N
            complete[row, col] = complete[col, row] = F(list[k += 1])
        end
    end
    complete
end

function Base.convert(::Type{Array{F,2}}, lm::PairwiseListMatrix{T, false, VT}) where {T,F,VT}
    N = lm.nelements
    complete = Array{F}(undef, N, N)
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

# Matrix(m::Symmetric{T, PairwiseListMatrix{T,D,VT}}) where {T,D,VT} = Matrix(m.data)

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
function PairwiseListMatrix(list::AbstractVector{T},
                           diagonal::Bool = false,
                           diagonalvalue::T = _diagonal_value(diagonal, T)) where T
    VT = typeof(list)
    if diagonal
        nelements = _nelements_with_diagonal(length(list))
        return PairwiseListMatrix{T, true, VT}(list, convert(VT,T[]), nelements)
    else
        nelements = _nelements(length(list))
        return PairwiseListMatrix{T, false, VT}(list,
            fill!(convert(VT, Array{T}(undef, nelements)), diagonalvalue),
            nelements)
    end
end

# AbstractArray methods
# =====================

Base.size(m::PairwiseListMatrix) = (m.nelements, m.nelements)
Base.length(m::PairwiseListMatrix) = m.nelements * m.nelements

Base.eltype(m::PairwiseListMatrix{T, diagonal, VT}) where {T, diagonal, VT} = T

for F in [:similar, :copy, :zeros, :ones]
    @eval begin
        function Base.$F(m::PairwiseListMatrix{T, diagonal, VT}) where {T, diagonal, VT}
            list = $F(m.list)
            VOUT = typeof(list)
            diag = convert(VOUT, $F(m.diag))
            PairwiseListMatrix{T, diagonal, VOUT}(list, diag, copy(m.nelements))
        end
    end
end

function Base.similar(m::PairwiseListMatrix{T, diagonal, VT},
                      ::Type{S}) where {T, diagonal, VT, S}
    list = similar(m.list, S)
    VOUT = typeof(list)
    diag = convert(VOUT, similar(m.diag, S))
    PairwiseListMatrix{S, diagonal, VOUT}(list, diag, copy(m.nelements))
end

# Indexing (getindex)
# ===================

function Base.getindex(lm::PairwiseListMatrix{T, true, VT}, i::Int, j::Int) where {T, VT}
    if i <= j
        return(lm.list[ij2k(i, j, lm.nelements, Val{true})])
    else
        return(lm.list[ij2k(j, i, lm.nelements, Val{true})])
    end
end

function Base.getindex(lm::PairwiseListMatrix{T, false, VT}, i::Int, j::Int) where {T, VT}
    if i < j
        return(lm.list[ij2k(i, j, lm.nelements, Val{false})])
    elseif i > j
        return(lm.list[ij2k(j, i, lm.nelements, Val{false})])
    else
        return(lm.diag[i])
    end
end

function Base.getindex(lm::PairwiseListMatrix{T, true, VT}, i::Int) where {T, VT}
    n = lm.nelements
    row = Int(ceil(i/n))
    col = i - n*(row-1)
    if row <= col
        @inbounds return lm.list[ij2k(row, col, n, Val{true})]
    else
        @inbounds return lm.list[ij2k(col, row, n, Val{true})]
    end
end

function Base.getindex(lm::PairwiseListMatrix{T, false, VT}, i::Int) where {T, VT}
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

function Base.setindex!(lm::PairwiseListMatrix{T, true, VT}, v, i::Int, j::Int) where {T, VT}
    if i <= j
        return setindex!(lm.list, v, ij2k(i, j, lm.nelements, Val{true}))
    else
        return setindex!(lm.list, v, ij2k(j, i, lm.nelements, Val{true}))
    end
end

function Base.setindex!(lm::PairwiseListMatrix{T, false, VT}, v, i::Int, j::Int) where {T, VT}
    if i < j
        return setindex!(lm.list, v, ij2k(i, j, lm.nelements, Val{false}))
    elseif i > j
        return setindex!(lm.list, v, ij2k(j, i, lm.nelements, Val{false}))
    else
        return setindex!(lm.diag, v, i)
    end
end

function Base.setindex!(lm::PairwiseListMatrix{T, true, VT}, v, i::Int) where {T, VT}
    n = lm.nelements
    row = Int(ceil(i/n))
    col = i - n*(row-1)
    if row <= col
        return setindex!(lm.list, v, ij2k(row, col, n, Val{true}))
    else
        return setindex!(lm.list, v, ij2k(col, row, n, Val{true}))
    end
end

function Base.setindex!(lm::PairwiseListMatrix{T, false, VT}, v, i::Int) where {T, VT}
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

LinearAlgebra.transpose(lm::PairwiseListMatrix)  = lm
LinearAlgebra.transpose!(lm::PairwiseListMatrix) = lm

LinearAlgebra.adjoint(lm::PairwiseListMatrix)  = lm
LinearAlgebra.adjoint!(lm::PairwiseListMatrix) = lm

# lu
# ==

function LinearAlgebra.lu!(lm::PairwiseListMatrix, args...; kargs...)
    LinearAlgebra.lu!(Matrix(lm), args...; kargs...)
end

# Diagonal
# ========

"""
Returns a vector of type `VT` from a `PairwiseListMatrix{T, false, VT}` that has
the diagonal values.

```jldoctest
julia> using PairwiseListMatrices

julia> plm = PairwiseListMatrix([10,20,30,40,50,60], true)
3×3 PairwiseListMatrix{Int64,true,Array{Int64,1}}:
 10  20  30
 20  40  50
 30  50  60

julia> diagonal(plm)
3-element Array{Int64,1}:
 10
 40
 60

```
"""
diagonal(lm::PairwiseListMatrix{T, false, VT}) where {T,VT} = lm.diag

function diagonal(lm::PairwiseListMatrix{T, true, VT}) where {T,VT}
    convert(VT, T[ lm[i,i] for i in 1:lm.nelements ])
end

# Customizing broadcasting
# ========================

function Base.BroadcastStyle(::Type{PairwiseListMatrix{T, D, VT}}) where {T, D, VT}
    Broadcast.DefaultArrayStyle{2}()
end

function Base.Broadcast.broadcasted(::Broadcast.DefaultArrayStyle{2},
                                    fun::F,
                                    plm::PairwiseListMatrix{T,D,V}) where {F<:Function,T,D,V}
    list = fun.(plm.list)
    E = eltype(list)
    if D
        diag = Array{E, 1}()
    else
        _diag = fun.(plm.diag)
        if eltype(_diag) !== E
            diag = copyto!(similar(_diag, E), _diag)
        else
            diag = _diag
        end
    end
    PairwiseListMatrix{E, D, typeof(list)}(list, diag, plm.nelements)
end

# Slow for multiple operations (it doesn't fuse them), e.g.: 0.5 .* $plm ./ 10.0
function Base.Broadcast.broadcasted(::Broadcast.DefaultArrayStyle{2},
                                    fun::F,
                                    plm::PairwiseListMatrix{T,D,V},
                                    n::Number) where {F<:Function,T,D,V}
    list = broadcast(fun, plm.list, n)
    E = eltype(list)
    if D
        diag = Array{E, 1}()
    else
        _diag = broadcast(fun, plm.diag, n)
        if eltype(_diag) !== E
            diag = copyto!(similar(_diag, E), _diag)
        else
            diag = _diag
        end
    end
    PairwiseListMatrix{E, D, typeof(list)}(list, diag, plm.nelements)
end

function Base.Broadcast.broadcasted(::Broadcast.DefaultArrayStyle{2},
                                    fun::F,
                                    n::Number,
                                    plm::PairwiseListMatrix{T,D,V}) where {F<:Function,T,D,V}
    list = broadcast(fun, n, plm.list)
    E = eltype(list)
    if D
        diag = Array{E, 1}()
    else
        _diag = broadcast(fun, n, plm.diag)
        if eltype(_diag) !== E
            diag = copyto!(similar(_diag, E), _diag)
        else
            diag = _diag
        end
    end
    PairwiseListMatrix{E, D, typeof(list)}(list, diag, plm.nelements)
end

@inline function Base.copyto!(dest::PairwiseListMatrix{T, D, VT},
                              bc::Base.Broadcast.Broadcasted{Nothing}) where {T, D, VT}
    axes(dest) == axes(bc) || Base.Broadcast.throwdm(axes(dest), axes(bc))
    bc_ = Base.Broadcast.preprocess(dest, bc)
    @iterateupper dest true list[k] = :($bc_)[CartesianIndex(i,j)] # slow: bc_ has a plm
    return dest
end

# Unary operations (faster than default methods)
# ==============================================

for F in (:abs, :-)
    @eval begin
        function Base.$F(lm::PairwiseListMatrix{T,true,VT}) where {T,VT}
            list = $F(lm.list)
            VOUT = typeof(list)
            diag = convert(VOUT, lm.diag)
            PairwiseListMatrix{eltype(list), true, VOUT}(list, diag, lm.nelements)
        end
        function Base.$F(lm::PairwiseListMatrix{T,false,VT}) where {T,VT}
            list = $F(lm.list)
            VOUT = typeof(list)
            diag = convert(VOUT, $F(lm.diag))
            PairwiseListMatrix{eltype(list), false, VOUT}(list, diag, lm.nelements)
        end
    end
end

for F in (:map, ) # , :broadcast)
    @eval begin
        function Base.$F(f, plm::PairwiseListMatrix{T,true,VT}) where {T,VT}
            list = $F(f, plm.list)
            VOUT = typeof(list)
            diag = convert(VOUT, plm.diag)
            PairwiseListMatrix{eltype(list), true, VOUT}(list, diag, plm.nelements)
        end

        function Base.$F(f, plm::PairwiseListMatrix{T,false,VT}) where {T,VT}
            list = $F(f, plm.list)
            VOUT = typeof(list)
            diag = convert(VOUT, $F(f, plm.diag))
            PairwiseListMatrix{eltype(list), false, VOUT}(list, diag, plm.nelements)
        end
    end
end

for F in (:map!, ) # , :broadcast!)
    @eval begin
        function Base.$F(f, plm::PairwiseListMatrix{T,true,VT}) where {T,VT}
            $F(f, plm.list)
            plm
        end

        function Base.$F(f, plm::PairwiseListMatrix{T,false,VT}) where {T,VT}
            $F(f, plm.list)
            $F(f, plm.diag)
            plm
        end
    end
end

LinearAlgebra.svd(m::PairwiseListMatrix; kargs...) = LinearAlgebra.svd(Matrix(m), kargs...)

function LinearAlgebra.svd(m::UpperTriangular{T,PairwiseListMatrix{T,D,VT}}) where {T,D,VT}
    svd(Matrix(m))
end

function LinearAlgebra.svd(m::LowerTriangular{T,PairwiseListMatrix{T,D,VT}}) where {T,D,VT}
    svd(Matrix(m))
end

# Faster mean and sum
# ===================

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

function _test_and_setsum(::Type{T}, dims, N) where T
    if dims == 1
        sum_i = zeros(T, 1, N)
    elseif dims == 2
        sum_i = zeros(T, N, 1)
    else
        throw(ErrorException("dims should be 1 or 2"))
    end
    sum_i
end

_sum(m::PairwiseListMatrix{T, false}) where {T} = T(2)*sum(m.list) + sum(m.diag)
_sum(m::PairwiseListMatrix{T, true}) where {T} = T(2)*sum(m.list) - sum(diagonal(m))

function Base.sum(lm::PairwiseListMatrix{T,true,VT}; dims::Union{Int,Colon}=:) where {T,VT}
    if isa(dims, Colon)
        return _sum(lm)
    end
    N = lm.nelements
    sum_i = _test_and_setsum(T, dims, N)
    _sum_kernel!(sum_i, lm.list, N)
end

function Base.sum(lm::PairwiseListMatrix{T,false,VT};
                  dims::Union{Int,Colon}=:) where {T, VT}
    if isa(dims, Colon)
        return _sum(lm)
    end
    N = lm.nelements
    sum_i = _test_and_setsum(T, dims, N)
    _sum_kernel!(sum_i, lm.diag, lm.list, N)
end

function Statistics.mean(m::PairwiseListMatrix; dims::Union{Int,Colon}=:)
    if isa(dims, Colon)
        sum(m)/length(m)
    else
        sum(m, dims=dims) ./ m.nelements
    end
end

# Sum/Mean without diagonal: sum/mean_nodiag
# ------------------------------------------

function _sum_nodiag_kernel!(sum_i, lm::PairwiseListMatrix{T,true,VT}, N) where {T,VT}
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

function _sum_nodiag_kernel!(sum_i, lm::PairwiseListMatrix{T,false,VT}, N) where {T,VT}
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

_sum_nodiag(m::PairwiseListMatrix{T,false,VT}) where {T,VT} = T(2) * sum(m.list)
function _sum_nodiag(m::PairwiseListMatrix{T,true,VT}) where {T,VT}
    T(2) * sum(m.list) - sum(diagonal(m))
end

"Sum the values outside the diagonal"
function sum_nodiag(lm::PairwiseListMatrix{T,diagonal,VT};
                    dims::Union{Int,Colon}=:) where {T,diagonal,VT}
    if isa(dims, Colon)
        return _sum_nodiag(lm)
    end
    N = lm.nelements
    sum_i = _test_and_setsum(T, dims, N)
    _sum_nodiag_kernel!(sum_i, lm, N)
end

"Mean of the values outside the diagonal"
function mean_nodiag(m::PairwiseListMatrix; dims::Union{Int,Colon}=:)
    if isa(dims, Colon)
        return sum_nodiag(m) / (length(m) - m.nelements)
    end
    sum_nodiag(m, dims=dims) ./ (m.nelements-1)
end

# Operations on Vector{PairwiseListMatrix}
# ========================================

# sum
# ---

@inline _has_diagonal(x::PairwiseListMatrix{T,diagonal,VT}) where {T,diagonal,VT} = diagonal

function Base.sum(list::Vector{PairwiseListMatrix{T,diagonal,VT}}) where {T,diagonal,VT}
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

# mean
# ----

function Statistics.mean(list::Vector{PairwiseListMatrix{T,diagonal,VT}}) where {T,diagonal,VT}
    sum(list) ./ length(list)
end

# std
# ---

function Statistics.varm(list::Vector{PairwiseListMatrix{T,D,VT}},
                   mean::PairwiseListMatrix{M,D,VM}) where {T,D,VT,M,VM}
    samples = length(list)
    if samples < 2
        throw(ErrorException("You need at least 2 samples."))
    end
    mean_list = mean.list
    N = length(mean_list)
    OutType = promote_type(T, M)
    out = fill!(PairwiseListMatrix{OutType,D,Vector{OutType}}(mean.nelements),zero(OutType))
    out_list = out.list
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

function Statistics.var(list::Vector{PairwiseListMatrix{T,D,VT}}; mean=nothing) where {T,D,VT}
    if mean === nothing
        varm(list, Statistics.mean(list))
    else
        varm(list, mean)
    end
end

function Statistics.std(list::Vector{PairwiseListMatrix{T,D,VT}}; mean=nothing) where {T,D,VT}
    sqrt.(var(list, mean=mean))
end

# zscore
# ======

"""
This function takes a vector of `PairwiseListMatrix` objects and a `PairwiseListMatrix` and
fill the matrix with the zscore value using the median and std of the vector.
"""
function zscore!(list::Vector{PairwiseListMatrix{L,D,VL}},
                 mat::PairwiseListMatrix{E,D,VE}) where {L,D,VL,E,VE}
    list_mean = mean(list)
    if size(mat) != size(list_mean)
        throw(ErrorException("PairwiseListMatrices must have the same size"))
    end
    list_std  = std(list, mean=list_mean)

    mat_list = mat.list
    std_list = list_std.list
    mean_list = list_mean.list

    @inbounds for i in 1:length(mat_list)
        s = std_list[i]
        m = mean_list[i]
        value = mat_list[i]
        if (m - value + one(E)) ≉ one(E)
            if (s + one(E)) ≉ one(E)
                mat_list[i] = (value - m)/s
            else
                mat_list[i] = NaN
            end
        else
            mat_list[i] = zero(E)
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
            if (m - value + one(E)) ≉ one(E)
                if (s + one(E)) ≉ one(E)
                    mat_diag[i] = (value - m)/s
                else
                    mat_diag[i] = NaN
                end
            else
                mat_diag[i] = zero(E)
            end
        end
    end

    mat
end

"It's like `zscore!` but without modifying the `PairwiseListMatrix`"
zscore(list, mat) = zscore!(list, copy(mat))

# NamedArrays
# ===========

# labels
# ------

"""
It gets the labels of a PairwiseListMatrix.

```jldoctest
julia> using PairwiseListMatrices

julia> plm  = PairwiseListMatrix(ones(3), false)
3×3 PairwiseListMatrix{Float64,false,Array{Float64,1}}:
 0.0  1.0  1.0
 1.0  0.0  1.0
 1.0  1.0  0.0

julia> getlabels(plm)
3-element Array{String,1}:
 "1"
 "2"
 "3"

julia> nplm  = setlabels(plm, ["a","b","c"])
3×3 Named PairwiseListMatrix{Float64,false,Array{Float64,1}}
A ╲ B │   a    b    c
──────┼──────────────
a     │ 0.0  1.0  1.0
b     │ 1.0  0.0  1.0
c     │ 1.0  1.0  0.0

julia> getlabels(nplm)
3-element Array{String,1}:
 "a"
 "b"
 "c"

```
"""
function getlabels(plm::PairwiseListMatrix{T,D,TV}) where {T,D,TV}
    String[ string(i) for i in 1:(plm.nelements) ]
end

getlabels(nplm::NamedArray{T,2,PairwiseListMatrix{T,D,TV},DN}) where {T,D,TV,DN} = names(nplm)[1]

"It changes the labels of a Named PairwiseListMatrix"
function setlabels!(nplm::NamedArray{T,2,PairwiseListMatrix{T,D,TV},DN},
                    labels::Vector{String}) where {T,D,TV,DN}
    NE = nplm.array.nelements
    NL = length(labels)
    @assert NL == NE "The number of elements, $NE, must be equal to the number of labels, $NL"
    setnames!(nplm, labels, 1)
    setnames!(nplm, labels, 2)
    nplm
end

"""
Creates a Named PairwiseListMatrix.

```jldoctest
julia> using PairwiseListMatrices

julia> plm  = PairwiseListMatrix(ones(3), false)
3×3 PairwiseListMatrix{Float64,false,Array{Float64,1}}:
 0.0  1.0  1.0
 1.0  0.0  1.0
 1.0  1.0  0.0

julia> nplm  = setlabels(plm, ["a","b","c"])
3×3 Named PairwiseListMatrix{Float64,false,Array{Float64,1}}
A ╲ B │   a    b    c
──────┼──────────────
a     │ 0.0  1.0  1.0
b     │ 1.0  0.0  1.0
c     │ 1.0  1.0  0.0

```
"""
function setlabels(nplm::NamedArray{T,2,PairwiseListMatrix{T,D,TV},DN},
                   labels::Vector{String}) where {T,D,TV,DN}
    setlabels!(copy(nplm), labels)
end

function setlabels(plm::PairwiseListMatrix{T,D,TV}, labels::Vector{String}) where {T,D,TV}
    nplm = NamedArray(plm)
    setlabels!(nplm, labels)
end

# Delegate to PairwiseListMatrix
# ------------------------------

for F in (:getlist, :getdiag, :lengthlist, :sum_nodiag, :mean_nodiag,
          :diagonal)
    @eval begin
        function $F(m::NamedArray{T,2,PairwiseListMatrix{T,D,TV},DN}) where {T,D,TV,DN}
            $F(m.array)
        end
    end
end

# Keep names
for F in (:map, :map!, :broadcast, :broadcast!)
    @eval begin
        function Base.$F(f, m::NamedArray{T,2,PairwiseListMatrix{T,D,TV},DN}) where {T,D,TV,DN}
            NamedArray(($F)(f, m.array), m.dicts, m.dimnames)
        end
    end
end

# Tables
# ======

"""
Creates a `Matrix{Any}`, labels are stored in the columns 1 and 2, and the
values in the column 3. Diagonal values are included by default.

```jldoctest
julia> using PairwiseListMatrices

julia> plm = PairwiseListMatrix([10,20,30], false)
3×3 PairwiseListMatrix{Int64,false,Array{Int64,1}}:
  0  10  20
 10   0  30
 20  30   0

julia> to_table(plm)
6×3 Array{Any,2}:
 "1"  "1"   0
 "1"  "2"  10
 "1"  "3"  20
 "2"  "2"   0
 "2"  "3"  30
 "3"  "3"   0

julia> to_table(plm, diagonal=false)
3×3 Array{Any,2}:
 "1"  "2"  10
 "1"  "3"  20
 "2"  "3"  30

```
"""
function to_table(plm::PairwiseListMatrix; diagonal::Bool = true, labels = getlabels(plm))
    N = plm.nelements
    table = Array{Any}(undef, diagonal ? div(N*(N+1),2) : div(N*(N-1),2), 3)
    t = 0
    @iterateupper plm diagonal begin
        :($t) += 1
        :($table)[:($t), 1] = :($labels)[i]
        :($table)[:($t), 2] = :($labels)[j]
        :($table)[:($t), 3] = list[k]
    end
    table
end

function to_table(nplm::NamedArray{T,2,PairwiseListMatrix{T,D,TV},DN};
                 diagonal::Bool = true,
                 labels::Vector{String} = getlabels(nplm)) where {T,D,TV,DN}
    to_table(nplm.array, diagonal=diagonal, labels=labels)
end


"""
It takes a `PairwiseListMatrix` and converts it to a `Dict` of `Symbol`s to arrays.
The returned dictionary can be easily converted into a `DataFrame`.

```jldoctest
julia> using PairwiseListMatrices, DataFrames

julia> nplm = setlabels(PairwiseListMatrix([10,20,30], false), ["a","b","c"])
3×3 Named PairwiseListMatrix{Int64,false,Array{Int64,1}}
A ╲ B │  a   b   c
──────┼───────────
a     │  0  10  20
b     │ 10   0  30
c     │ 20  30   0

julia> dict = to_dict(nplm, diagonal=false)
Dict{Symbol,Array{T,1} where T} with 3 entries:
  :values => [10, 20, 30]
  :j      => ["b", "c", "c"]
  :i      => ["a", "a", "b"]

julia> DataFrame(dict)
3×3 DataFrames.DataFrame
│ Row │ i      │ j      │ values │
│     │ String │ String │ Int64  │
├─────┼────────┼────────┼────────┤
│ 1   │ a      │ b      │ 10     │
│ 2   │ a      │ c      │ 20     │
│ 3   │ b      │ c      │ 30     │
```
"""
function to_dict(plm::PairwiseListMatrix{T,D,TV};
                 diagonal::Bool=true,
                 labels::Vector{String} = getlabels(plm)) where {T,D,TV}
    N = plm.nelements
    L = diagonal ? div(N*(N+1),2) : div(N*(N-1),2)
    I = Array{String}(undef, L)
    J = Array{String}(undef, L)
    K = Array{T}(undef, L)
    t = 0
    @iterateupper plm diagonal begin
        :($t) += 1
        :($I)[:($t)] = :($labels)[i]
        :($J)[:($t)] = :($labels)[j]
        :($K)[:($t)] = list[k]
    end
    Dict(:i => I, :j => J, :values => K)
end

function to_dict(nplm::NamedArray{T,2,PairwiseListMatrix{T,D,TV},DN};
                 diagonal::Bool = true,
                 labels::Vector{String} = getlabels(nplm)) where {T,D,TV,DN}
    to_dict(nplm.array, diagonal=diagonal, labels=labels)
end

"""
Creation of a `PairwiseListMatrix` from a `Matrix`, `DataFrame` or similar structure.
By default the columns with the labels for i (slow) and j (fast) are 1 and 2.
Values are taken from the column 3 by default.

```jldoctest
julia> using PairwiseListMatrices, Pkg, DelimitedFiles

julia> import PairwiseListMatrices

julia> filename = joinpath(dirname(pathof(PairwiseListMatrices)), "..", "test", "example.csv");

julia> dat = readdlm(filename, ',')
3×3 Array{Any,2}:
 "A"  "B"  10
 "A"  "C"  20
 "B"  "C"  30

julia> from_table(dat, false)
3×3 Named PairwiseListMatrix{Any,false,Array{Any,1}}
A ╲ B │       A        B        C
──────┼──────────────────────────
A     │ nothing       10       20
B     │      10  nothing       30
C     │      20       30  nothing
```

This is also useful to create a `PairwiseListMatrix` from a `DataFrame`:

```
julia> using PairwiseListMatrices, DataFrames, CSV, Pkg

julia> import PairwiseListMatrices

julia> filename = joinpath(dirname(pathof(PairwiseListMatrices)), "..", "test", "example.csv");

julia> df = CSV.read(filename, header=false)
3×3 DataFrames.DataFrame
│ Row │ Column1 │ Column2 │ Column3 │
│     │ String⍰ │ String⍰ │ Int64⍰  │
├─────┼─────────┼─────────┼─────────┤
│ 1   │ A       │ B       │ 10      │
│ 2   │ A       │ C       │ 20      │
│ 3   │ B       │ C       │ 30      │

julia> from_table(df, false)
3×3 Named PairwiseListMatrix{Union{Missing, Int64},false,Array{Union{Missing, Int64},1}}
A ╲ B │       A        B        C
──────┼──────────────────────────
A     │ missing       10       20
B     │      10  missing       30
C     │      20       30  missing
```
"""
function from_table(table,
                    diagonal::Bool;
                    labelcols::Vector{Int} = [1,2],
                    valuecol::Int = 3,
                    diagonalvalue = :default)
    @assert size(table,2) >= 3
    values = table[:,valuecol]
    if diagonalvalue == :default
        diagonalvalue = _diagonal_value(diagonal, eltype(values))
    end
    plm  = PairwiseListMatrix(values, diagonal, diagonalvalue)
    nplm = NamedArray(plm)
    if length(labelcols) == 2
        unique_labels = unique(vcat(table[:, labelcols[1]],
                                    table[:, labelcols[2]]))
        setlabels!(nplm, String[ string(lab) for lab in unique_labels ])
    end
    nplm
end

# write...
# ========

"""
This function takes the filename as first argument and a `PairwiseListMatrix` as second
argument. If the `diagonal` keyword argument is `true` (default), the diagonal is included
in the output. The keyword argument `delim` (by default is `'\t'`) allows to modified the
character used as delimiter.

```jldoctest
julia> using PairwiseListMatrices, DelimitedFiles

julia> plm  = PairwiseListMatrix(trues(3), false)
3×3 PairwiseListMatrix{Bool,false,BitArray{1}}:
 false   true   true
  true  false   true
  true   true  false

julia> writedlm("example.csv", plm, diagonal=false, delim=',')

julia> println(read("example.csv", String))
1,2,true
1,3,true
2,3,true

```
"""
function DelimitedFiles.writedlm(filename::String,
                       plm::PairwiseListMatrix{T,D,TV};
                       diagonal::Bool = true,
                       delim::Char = '\t',
                       labels::Vector{String} = getlabels(plm)) where {T,D,TV}
    open(filename, "w") do fh
        @iterateupper plm diagonal println(:($fh), :($labels)[i], :($delim), :($labels)[j], :($delim), list[k])
    end
end

function DelimitedFiles.writedlm(filename::String,
                       nplm::NamedArray{T,2,PairwiseListMatrix{T,D,TV},DN};
                       diagonal::Bool = true,
                       delim::Char = '\t',
                       labels::Vector{String} = getlabels(nplm)) where {T,D,TV,DN}
    writedlm(filename, nplm.array, diagonal=diagonal, delim=delim, labels=labels)
end

# triu! & triu
# ============

function LinearAlgebra.triu!(mat::PairwiseListMatrix)
    throw(ErrorException("PairwiseListMatrix must be Symmetric, use triu instead of triu!"))
end
LinearAlgebra.triu!(mat::PairwiseListMatrix, k::Int) = triu!(mat)
LinearAlgebra.triu(mat::PairwiseListMatrix) = triu(Matrix(mat))
LinearAlgebra.triu(mat::PairwiseListMatrix, k::Int) = triu(Matrix(mat), k::Int)

# join
# ====

"""
This function join two PairwiseListMatrices by their labels, returning two
PairwiseListMatrices with same size and labels. There are 4 `kind`s of joins:

- `:inner` : Intersect. The output matrices only include the labels that are in both PairwiseListMatrices
- `:outer` : Union. Include the labels of the two PairwiseListMatrices.
- `:left` : Only use labels from the first argument.
- `:right` : Only use labels from the second argument.

`NaN`s are filled in where needed to complete joins. The default value for missing values
can be changed passing a tuple to `missing`.

```jldoctest
julia> using PairwiseListMatrices

julia> l = setlabels(PairwiseListMatrix([1.,2.,3.], false), ["a","b","c"]) # a b c
3×3 Named PairwiseListMatrix{Float64,false,Array{Float64,1}}
A ╲ B │   a    b    c
──────┼──────────────
a     │ 0.0  1.0  2.0
b     │ 1.0  0.0  3.0
c     │ 2.0  3.0  0.0

julia> r = setlabels(PairwiseListMatrix([1.,2.,3.], false), ["b","c","d"]) # b c d
3×3 Named PairwiseListMatrix{Float64,false,Array{Float64,1}}
A ╲ B │   b    c    d
──────┼──────────────
b     │ 0.0  1.0  2.0
c     │ 1.0  0.0  3.0
d     │ 2.0  3.0  0.0

julia> join(l, r, kind=:inner) # b c
(2×2 Named PairwiseListMatrix{Float64,false,Array{Float64,1}}
A ╲ B │   b    c
──────┼─────────
b     │ 0.0  3.0
c     │ 3.0  0.0, 2×2 Named PairwiseListMatrix{Float64,false,Array{Float64,1}}
A ╲ B │   b    c
──────┼─────────
b     │ 0.0  1.0
c     │ 1.0  0.0)

julia> join(l, r, kind=:outer) # a b c d
(4×4 Named PairwiseListMatrix{Float64,false,Array{Float64,1}}
A ╲ B │   a    b    c    d
──────┼───────────────────
a     │ 0.0  1.0  2.0  NaN
b     │ 1.0  0.0  3.0  NaN
c     │ 2.0  3.0  0.0  NaN
d     │ NaN  NaN  NaN  NaN, 4×4 Named PairwiseListMatrix{Float64,false,Array{Float64,1}}
A ╲ B │   a    b    c    d
──────┼───────────────────
a     │ NaN  NaN  NaN  NaN
b     │ NaN  0.0  1.0  2.0
c     │ NaN  1.0  0.0  3.0
d     │ NaN  2.0  3.0  0.0)

```
"""
function Base.join(
        left::NamedArray{L,2,PairwiseListMatrix{L,DL,VL},NL},
        right::NamedArray{R,2,PairwiseListMatrix{L,DR,VR},NR};
        kind::Symbol = :inner,
        missing::Tuple{L,R} = (L(NaN), R(NaN))
        ) where {L <: AbstractFloat, R <: AbstractFloat,DL,DR,VL,VR,NL,NR}

    labels_left = getlabels(left)
    labels_right = getlabels(right)

    if labels_left == labels_right
        return(left, right)
    end

    if kind == :inner
        out_labels = intersect(labels_left, labels_right)
        N = length(out_labels)
        out_L = PairwiseListMatrix(L, N, DL)
        out_R = PairwiseListMatrix(R, N, DR)
        for i in 1:N
            li = out_labels[i]
            for j in i:N
                lj = out_labels[j]
                out_L[i,j] = left[li, lj]
                out_R[i,j] = right[li, lj]
            end
        end
        return(setlabels(out_L, out_labels), setlabels(out_R, out_labels))
    elseif kind == :left
        out_labels = labels_left
        N = length(out_labels)
        out_R = PairwiseListMatrix(R, N, DR)
        for i in 1:N
            li = out_labels[i]
            flag_i_r = li in labels_right
            for j in i:N
                lj = out_labels[j]
                out_R[i,j] = (flag_i_r && (lj in labels_right)) ? right[li, lj] : missing[2]
            end
        end
        return(left, setlabels(out_R, out_labels))
    elseif kind == :right
        out_labels = labels_right
        N = length(out_labels)
        out_L = PairwiseListMatrix(L, N, DL)
        for i in 1:N
            li = out_labels[i]
            flag_i_l = li in labels_left
            for j in i:N
                lj = out_labels[j]
                out_L[i,j] = (flag_i_l  && (lj in labels_left)) ? left[li, lj] : missing[1]
            end
        end
        return(setlabels(out_L, out_labels), right)
    elseif kind == :outer
        out_labels = union(labels_left, labels_right)
        N = length(out_labels)
        out_L = PairwiseListMatrix(L, N, DL)
        out_R = PairwiseListMatrix(R, N, DR)
        for i in 1:N
            li = out_labels[i]
            flag_i_l = li in labels_left
            flag_i_r = li in labels_right
            for j in i:N
                lj = out_labels[j]
                out_L[i,j] = (flag_i_l && (lj in labels_left))  ? left[li, lj] : missing[1]
                out_R[i,j] = (flag_i_r && (lj in labels_right)) ? right[li, lj] : missing[2]
            end
        end
        return(setlabels(out_L , out_labels), setlabels(out_R, out_labels))
    else
        throw(ArgumentError("Unknown kind of join requested: use :inner, :left, :right or :outer"))
    end

end
