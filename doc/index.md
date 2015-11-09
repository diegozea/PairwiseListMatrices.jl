---
order: 0
title: PairwiseListMatrices
author: Diego Javier Zea
modules: [PairwiseListMatrices]
...

# Installation

This package is registered on Julia's METADATA.jl and it can be installed through the Julia's REPL:
```{.julia execute="false"}
Pkg.add("PairwiseListMatrices")
```
If the package is installed on your system, you can load it with:
```julia
using PairwiseListMatrices
```

# Introduction

**PairwiseListMatrices** allows you to represent a (squared) **symmetric matrix** as a list of the values in the upper or lower triangular part of the matrix.
This matrices are common for representing **pairwise measures/comparisons** between the elements of one group when the used metric/distance satisfies the symmetry condition.
Also the **adjacency matrices of undirected graphs** can be represented with this kind of list/matrices.

## Simple example
The following symmetric matrix has 9 values. Their values could be thinked as pairwise measures between 3 elements:
```{.julia execute="false"}
julia> matrix = [  0 10 20
                  10  0 30
                  20 30  0 ]
3x3 Array{Int64,2}:
  0  10  20
 10   0  30
 20  30   0
```

Since all the diagonal members are zeros, this matrix could be represented as a vector/list of the 3 values on the triangular part:
```{.julia execute="false"}
julia> list = [10, 20, 30]
3-element Array{Int64,1}:
 10
 20
 30
```

The type `PairwiseListMatrix`, defined on this module, can be used for working with the list as a full symmetric matrix.
```{.julia execute="false"}
julia> plm = PairwiseListMatrix(list)
3x3 PairwiseListMatrices.PairwiseListMatrix{Int64,false}:
  0  10  20
 10   0  30
 20  30   0
```

# Implementation
If you are performing pairwise measures over `N` elements, storing all the `N*N` values on a `Matrix{T}` represents `sizeof(T)*(N*N)` bytes on memory. However, the lower and upper triangular part of the matrix are identical and could be stored in a single list. In this way, you are storing the green value only once:

<div style="text-align: center"><img width=300px src="PLM.svg"></div>

The diagonal values should be stored, since they could change at any time (i.e. yellow value). So you need `sizeof(T)*(N)` bytes for storing the diagonal values on a vector and `sizeof(T)*(N*(N-1))/2` bytes for storing the lower or upper triangular part of the matrix.
The type `PairwiseListMatrix{T, diagonal}` represents the symmetric matrix using only `sizeof(T)*(N*(N+1))/2` bytes instead of `sizeof(T)*(N*N)` bytes, saving almost 50% of the memory (the percent depends on `N`):

```julia
using Gadfly
plot(N -> 100.0 - ( 100.0 * div(N*(N+1), 2) / (N*N) ), 2, 500, Guide.xlabel("N"), Guide.ylabel("% of saved memory"))
```
As you can see in the schematic diagram, the difference between `PairwiseListMatrix{T, true}` and `PairwiseListMatrix{T, false}` is where the diagonal values are stored.
All `PairwiseListMatrix{T, diagonal}` have a `list` field for storing the values. If `diagonal` is true, the diagonal values are included on the `list` (i.e. yellow value) and the `diag` is empty. But if the `diagonal` value is `false` the diagonal values are stored on the `diag` vector.

```{.julia execute="false"}
type PairwiseListMatrix{T, diagonal} <: AbstractArray{T, 2}
  list::Vector{T}
  diag::Vector{T}
  labels::IndexedArray
  nelements::Int
end
```
The number of elements in the pairwise measure/comparisons or the number of nodes in the undirected graph is stored in `nelements` and used on indexing operations. This allows you to index the object like any other matrix.
The field `labels` is an Indexed Array, allowing to easily index the matrix using the labels.
For example, using the matrix of the figure and storing the diagonal values on the list:

```{.julia execute="false"}
julia> using PairwiseListMatrices

julia> plm = PairwiseListMatrix([1,1,0,1,0,1,1,0,0,0], ['A','B','C','D'], true)
4x4 PairwiseListMatrices.PairwiseListMatrix{Int64,true}:
 1  1  0  1
 1  0  1  1
 0  1  0  0
 1  1  0  0

julia> dump(plm)
PairwiseListMatrices.PairwiseListMatrix{Int64,true}
  list: Array(Int64,(10,)) [1,1,0,1,0,1,1,0,0,0]
  diag: Array(Int64,(0,)) Int64[]
  labels: IndexedArrays.IndexedArray{Char}
    items: Array(Char,(4,)) ['A','B','C','D']
    lookup: Dict{Char,Int64} len 4
      D: Int64 4
      B: Int64 2
      C: Int64 3
      A: Int64 1
  nelements: Int64 4
```

You can also create the matrix with the list without the diagonal values and fill the diagonal values after that:

```{.julia execute="false"}
julia> plm = PairwiseListMatrix([1,0,1,1,1,0], ['A','B','C','D'], false)
4x4 PairwiseListMatrices.PairwiseListMatrix{Int64,false}:
 0  1  0  1
 1  0  1  1
 0  1  0  0
 1  1  0  0

julia> setlabel!(plm, 1, 'A', 'A') # Set the yellow value using the labels
4-element Array{Int64,1}:
 1
 0
 0
 0

julia> plm
4x4 PairwiseListMatrices.PairwiseListMatrix{Int64,false}:
 1  1  0  1
 1  0  1  1
 0  1  0  0
 1  1  0  0

julia> dump(plm)
PairwiseListMatrices.PairwiseListMatrix{Int64,false}
  list: Array(Int64,(6,)) [1,0,1,1,1,0]
  diag: Array(Int64,(4,)) [1,0,0,0]
  labels: IndexedArrays.IndexedArray{Char}
    items: Array(Char,(4,)) ['A','B','C','D']
    lookup: Dict{Char,Int64} len 4
      D: Int64 4
      B: Int64 2
      C: Int64 3
      A: Int64 1
  nelements: Int64 4

```
# Ploting

The function `protovis` provides an **arc diagram** (since this could be a representation for an **adjacency matrix/list** of an undirected graph) and a **matrix visualization** on the web browser using [Protovis](http://mbostock.github.io/protovis/).

<div style="text-align: center"><img width=650px src="https://raw.githubusercontent.com/diegozea/PairwiseListMatrices.jl/master/doc/protovis_example.png"></div>

# Benchmark

`PairwiseListMatrix` is faster than a full matrix to make operation like `sum` and `mean` in the whole matrix, since it is cache efficient. However it is slower than a full matrix for reducing along dimensions.

 - [Creation benchmark](http://nbviewer.ipython.org/github/diegozea/PairwiseListMatrices.jl/blob/master/test/creation_bech.ipynb)
 - [Statistics benchmark](http://nbviewer.ipython.org/github/diegozea/PairwiseListMatrices.jl/blob/master/test/stats_bench.ipynb)
