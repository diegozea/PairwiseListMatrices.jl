"""
The `apply2list` function applies the function `fun`, the first argument, for 
each element on `list` but avoiding `getfield` calls inside the loop. 
The second argument of the macro is the `PairwiseListMatrix` that will be 
iterated. `fun` should take two arguments; the first is the `list` field of 
the `PairwiseListMatrix`, and the second is the index over that list at a 
given iteration.

```jldoctest
julia> using PairwiseListMatrices

julia> plm = PairwiseListMatrix([1,2,3], false)
3×3 PairwiseListMatrix{Int64,false,Array{Int64,1}}:
 0  1  2
 1  0  3
 2  3  0

julia> apply2list((list, k) -> println(list[k]), plm)
1
2
3

```
"""
function apply2list(fun, plm)
    list = plm.list
    for k in 1:length(list)
        fun(list, k)
    end
end

"""
The `apply2diag` function applies the function `fun`, the first argument, 
for each element on the `diag` field of a`PairwiseListMatrix{T,false,VT}` 
but avoiding `getfield` calls inside the loop. The second argument of the 
macro is the `PairwiseListMatrix` that will be iterated. `fun` should take 
two arguments; the first is the `diag` field of the `PairwiseListMatrix`, 
and the second is the index over that vector of diagonal elements at a 
given iteration.

```jldoctest
julia> using PairwiseListMatrices

julia> plm = PairwiseListMatrix([1,2,3], false)
3×3 PairwiseListMatrix{Int64,false,Array{Int64,1}}:
 0  1  2
 1  0  3
 2  3  0

julia>apply2diag((diag, k) -> diag[k] += 10k, plm)

julia> plm
3×3 PairwiseListMatrix{Int64,false,Array{Int64,1}}:
 10   1   2
  1  20   3
  2   3  30

```
"""
function apply2diag(fun, plm::PairwiseListMatrix{T,false,VT}) where {T, VT}
    diag = plm.diag
    for k in 1:length(diag)
        fun(diag, k)
    end
end

"""
The `apply2upper` function applies the `fun` function over the upper triangular 
part of the `PairwiseListMatrix`. Set `use_diag` to `true` if the diagonal 
needs to be included in the iteration. The function should take four 
arguments, first the `list` and diag fields of the `PairwiseListMatrix`, 
the second is the index `k` over that `list`, and the third and fourth are 
the `i` and `j` indexes for the position `k` in the upper triangular 
part of the matrix. 

```jldoctest
julia> using PairwiseListMatrices

julia> plm = PairwiseListMatrix([1,2,3], true)
2×2 PairwiseListMatrix{Int64,true,Array{Int64,1}}:
 1  2
 2  3

julia> mat = zeros(Int, 2, 2)
2×2 Array{Int64,2}:
 0  0
 0  0

julia> apply2upper((list, k, i, j) -> mat[i, j] = list[k], plm; use_diag = true)

julia> mat
2×2 Array{Int64,2}:
 1  2
 0  3

```
"""
function apply2upper(fun, plm::PairwiseListMatrix{T,false,VT}; use_diag::Bool=false) where {T, VT}
    N = plm.nelements
    if use_diag
        k = 0
        diag = plm.diag
        list = plm.list
        for i in 1:N
            for j in i:N
                if i != j
                    k += 1
                    fun(list, k, i, j)
                else
                    let list = diag, k = i
                        fun(list, k, i, j)
                    end
                end
            end
        end
    else
        k = 0
        list = plm.list
        for i in 1:(N-1)
            for j in (i+1):N
                k += 1
                fun(list, k, i, j)
            end
        end
    end
end

function apply2upper(fun, plm::PairwiseListMatrix{T,true,VT}; use_diag::Bool=false) where {T, VT}
    N = plm.nelements
    if hasdiagonal(plm)
        if use_diag
            k = 0
            list = plm.list
            for i in 1:N
                for j in i:N
                    k += 1
                    fun(list, k, i, j)
                end
            end
        else
            k = 0
            list = plm.list
            for i in 1:N
                for j in i:N
                    k += 1
                    if i != j
                        fun(list, k, i, j)
                    end
                end
            end
        end
    end
end
