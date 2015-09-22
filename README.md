# PairwiseListMatrices

Linux, OSX: [![Build Status](https://travis-ci.org/diegozea/PairwiseListMatrices.jl.svg?branch=master)](https://travis-ci.org/diegozea/PairwiseListMatrices.jl)

Windows: [![Build status](https://ci.appveyor.com/api/projects/status/p96sso5b23gi85mg/branch/master?svg=true)](https://ci.appveyor.com/project/diegozea/pairwiselistmatrices-jl/branch/master)

Code Coverage: [![Coverage Status](https://coveralls.io/repos/diegozea/PairwiseListMatrices.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/diegozea/PairwiseListMatrices.jl?branch=master) [![codecov.io](http://codecov.io/github/diegozea/PairwiseListMatrices.jl/coverage.svg?branch=master)](http://codecov.io/github/diegozea/PairwiseListMatrices.jl?branch=master)

This package allows you to use a paired list as a Matrix.  

`PairwiseListMatrix{T, L, diagonal}` is a (squared) symmetric matrix that stores a `list` of values of type `T` for the pairwise comparison/evaluation of `nelements`.
If `diagonal` is `true` the first element of the list is `1, 1` otherwise is `1, 2`.
If `diagonal` is `false`, the diagonal values are stored in a vector on the `diag` field.
Labels can be stored on the field `labels` as an `IndexedArray`.  

In pairwise calculations like `cor()` results are saved as `PairwiseListMatrix`, that is `N(N-1)/2` in space instead of `N*N`. This is useful to compare a large number of elements.

`PairwiseListMatrix` gives the option of save labels and allows to use them for indexing.

## Example

```julia
shell> cat example.csv
A,B,10
A,C,20
B,C,30

julia> data = readcsv("example.csv")
3x3 Array{Any,2}:
 "A"  "B"  10
 "A"  "C"  20
 "B"  "C"  30

julia> list = from_table(data, Int, ASCIIString, false)
3x3 PairwiseListMatrices.PairwiseListMatrix{Int64,ASCIIString,false}:
  0  10  20
 10   0  30
 20  30   0

julia> list.list
3-element Array{Int64,1}:
 10
 20
 30

julia> labels(list)
3-element IndexedArrays.IndexedArray{ASCIIString}:
 "A"
 "B"
 "C"

julia> getlabel(list, "A", "B")
10

julia> list[1,2]
10

julia> full(list)
3x3 Array{Int64,2}:
  0  10  20
 10   0  30
 20  30   0

```

## Benchmark
`PairwiseListMatrix` is faster than a full matrix to make operatation like `sum` and `mean` in the whole matrix, since is cache efficient. However it is slower than a full matrix for reducing along dimensions.

 - [Creation benchmark](http://nbviewer.ipython.org/github/diegozea/PairwiseListMatrices.jl/blob/master/test/creation_bech.ipynb)
 - [Statistics benchmark](http://nbviewer.ipython.org/github/diegozea/PairwiseListMatrices.jl/blob/master/test/stats_bench.ipynb)
