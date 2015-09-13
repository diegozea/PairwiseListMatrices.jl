# PairedListMatrices

Linux, OSX: [![Build Status](https://travis-ci.org/diegozea/PairedListMatrices.jl.svg?branch=master)](https://travis-ci.org/diegozea/PairedListMatrices.jl)

Windows: [![Build status](https://ci.appveyor.com/api/projects/status/u8ayy5ep2tnncnyp/branch/master?svg=true)](https://ci.appveyor.com/project/diegozea/pairedlistmatrices-jl/branch/master)

Code Coverage: [![Coverage Status](https://coveralls.io/repos/diegozea/PairedListMatrices.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/diegozea/PairedListMatrices.jl?branch=master) [![codecov.io](http://codecov.io/github/diegozea/PairedListMatrices.jl/coverage.svg?branch=master)](http://codecov.io/github/diegozea/PairedListMatrices.jl?branch=master)

This package allows you to use a paired list as a Matrix.
In pairwise calculations like `cor()`, saving the result as a `PairedListSymmetric` is `N(N-1)/2` in space, instead of `N*N`. This is very useful when you need to compare a large number of vectors:

```julia
julia> n = 60_000;

julia> Array(Float64, n, n);
ERROR: OutOfMemoryError()
 in call at essentials.jl:202
 in eval_user_input at REPL.jl:64
 [inlined code] from REPL.jl:93
 in anonymous at task.jl:68

julia> Array(Float64, div(n*(n-1),2));

julia>

```

`PairedListMatrices` is faster than other alternatives, since is cache efficient.
Also it gives you the option of save your labels and allows you to use them for indexing.

## Example

```julia
julia> using PairedListMatrices

julia> points = [ rand(3) for n in 1:3 ]
3-element Array{Array{Float64,1},1}:
 [0.12202515834618777,0.33975381354416134,0.7899163980736033]
 [0.0703637842317335,0.008293842156420261,0.8338452252525557]
 [0.18788223427528283,0.5878742934201919,0.23343241338325327]

julia> labels = ["O", "P", "Q"]
3-element Array{ASCIIString,1}:
 "O"
 "P"
 "Q"

julia> list = Array(Float64, 3);

julia> k=0;

julia> for i in 1:2
       for j in (i+1):3
       k += 1
       list[k] = cor(points[i], points[j])
       println(labels[i], " ", labels[j], " ", list[k])
       end
       end
O P 0.923814951778981
O Q -0.09394903823267092
P Q -0.46793753775974317

julia> mat = PairedListSymmetric(list, labels, 1.0)
3x3 PairedListMatrices.PairedListSymmetric{Float64,ASCIIString}:
  1.0        0.923815  -0.093949
  0.923815   1.0       -0.467938
 -0.093949  -0.467938   1.0

julia> mat.list
3-element Array{Float64,1}:
  0.923815
 -0.093949
 -0.467938

julia> mat.labels
3-element IndexedArrays.IndexedArray{ASCIIString}:
 "O"
 "P"
 "Q"

julia> getlabel(mat, "P", "Q")
-0.46793753775974317

```

## Benchmark

A `PairedListDiagonalSymmetric` (is diagonal because the list also includes the diagonal values) used as a matrix (`pairedlist_matindex`) can be filled faster the a Symmetric full matrix (`using_full_symmetric`). Is also faster than filling the full matrix as in `pairwise!` of **Distances.jl** (`distances_pairwise`). Filling a list (vector) and create a `PairedListDiagonalSymmetric` with it is also fast (`pairedlist_listindex`). When space is a problem `PairedListDiagonalSymmetric` everything is faster than saving the values on a sparse matrix (`using_full_symmetric`). The bechmark code is in the test folder.

```julia
julia> function pairedlist_matindex(vecs)
         n = length(vecs)
         list = PairedListDiagonalSymmetric(Float64, n)
         @inbounds for i in 1:n
           vec_i = vecs[i]
           for j in i:n
             list[i,j] = cor(vec_i, vecs[j])
           end
         end
         list
       end
pairedlist_matindex (generic function with 1 method)

julia> function pairedlist_listindex(vecs)
         n = length(vecs)
         list = Array(Float64, div(n*(n+1),2))
         k = 1
         @inbounds for i in 1:n
           vec_i = vecs[i]
           for j in i:n
             list[k] = cor(vec_i, vecs[j])
             k += 1
           end
         end
         PairedListDiagonalSymmetric(list)
       end
pairedlist_listindex (generic function with 1 method)

```

#### Small test set: [ rand(3) for n in 1:100 ]


**pairedlist_matindex**
Time per evaluation: 1.29 ms [870.25 μs, 1.70 ms]

using_full_symmetric
Time per evaluation: 1.43 ms [1.05 ms, 1.80 ms]

**pairedlist_listindex**
Time per evaluation: 1.56 ms [1.17 ms, 1.95 ms]

distances_pairwise
Time per evaluation: 1.61 ms [1.07 ms, 2.14 ms]

using_sparse_symmetric
Time per evaluation: 3.75 ms [3.35 ms, 4.14 ms]

#### Big test set: [ rand(3) for n in 1:1000 ]

**pairedlist_matindex**
Time per evaluation: 114.25 ms [105.46 ms, 123.04 ms]

**pairedlist_listindex**
Time per evaluation: 149.55 ms [139.16 ms, 159.94 ms]

using_full_symmetric
Time per evaluation: 154.07 ms [140.74 ms, 167.40 ms]

distances_pairwise
Time per evaluation: 161.17 ms [145.42 ms, 176.91 ms]

using_sparse_symmetric
Time per evaluation: 26.46 s
