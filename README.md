# PairwiseListMatrices

Julia 0.4: [![PairwiseListMatrices](http://pkg.julialang.org/badges/PairwiseListMatrices_0.4.svg)](http://pkg.julialang.org/?pkg=PairwiseListMatrices)
Julia 0.5: [![PairwiseListMatrices](http://pkg.julialang.org/badges/PairwiseListMatrices_0.5.svg)](http://pkg.julialang.org/?pkg=PairwiseListMatrices)  
Julia 0.6: [![PairwiseListMatrices](http://pkg.julialang.org/badges/PairwiseListMatrices_0.6.svg)](http://pkg.julialang.org/?pkg=PairwiseListMatrices)  

Linux, OSX: [![Build Status](https://travis-ci.org/diegozea/PairwiseListMatrices.jl.svg?branch=master)](https://travis-ci.org/diegozea/PairwiseListMatrices.jl)  
Windows: [![Build status](https://ci.appveyor.com/api/projects/status/p96sso5b23gi85mg/branch/master?svg=true)](https://ci.appveyor.com/project/diegozea/pairwiselistmatrices-jl/branch/master)  

Code Coverage: [![Coverage Status](https://coveralls.io/repos/diegozea/PairwiseListMatrices.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/diegozea/PairwiseListMatrices.jl?branch=master) [![codecov.io](http://codecov.io/github/diegozea/PairwiseListMatrices.jl/coverage.svg?branch=master)](http://codecov.io/github/diegozea/PairwiseListMatrices.jl?branch=master)

## Documentation  

[![stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://diegozea.github.io/PairwiseListMatrices.jl/stable)  
[![latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://diegozea.github.io/PairwiseListMatrices.jl/latest)  

## Description

This package allows you to use a pairwise list as a matrix:

![PLM](https://raw.githubusercontent.com/diegozea/PairwiseListMatrices.jl/master/docs/src/PLM_README.png)

```julia
mutable struct PairwiseListMatrix{T,diagonal,VT} <: AbstractArray{T, 2}
    list::VT
    diag::VT
    nelements::Int
    ...
end
```   

`PairwiseListMatrix{T, diagonal, VT}` is a (squared) symmetric matrix that stores a `list`
of type `VT` with values of type `T` for the pairwise comparison/evaluation of `nelements`.
If `diagonal` is `true` the first element of the list is `1, 1` otherwise is `1, 2`.
If `diagonal` is `false` the diagonal values are stored in a vector on the `diag` field.  

## Features  

#### Space  

In pairwise calculations like `cor()` if results are saved as `PairwiseListMatrix` the
space is `N(N+1)/2` instead of `N*N`. This is useful to compare a large number of elements,
because you are **saving ~ 50% of the memory.**  

#### Time  

`PairwiseListMatrix` is **faster than a full matrix** to make operatation like `sum` and
`mean` in the whole matrix, since it is cache efficient. However it is slower than a full
matrix for reducing along dimensions.  

 - [Creation benchmark](http://nbviewer.ipython.org/github/diegozea/PairwiseListMatrices.jl/blob/master/test/creation_bech.ipynb)
 - [Statistics benchmark](http://nbviewer.ipython.org/github/diegozea/PairwiseListMatrices.jl/blob/master/test/stats_bench.ipynb)

## Example

```julia
julia> # Pkg.add("PairwiseListMatrices")

julia> using PairwiseListMatrices

julia> plm  = PairwiseListMatrix([1,2,3], false)
3×3 PairwiseListMatrices.PairwiseListMatrix{Int64,false,Array{Int64,1}}:
 0  1  2
 1  0  3
 2  3  0

julia> nplm  = setlabels(plm, ["a","b","c"])
3×3 Named PairwiseListMatrices.PairwiseListMatrix{Int64,false,Array{Int64,1}}
A ╲ B │ a  b  c
──────┼────────
a     │ 0  1  2
b     │ 1  0  3
c     │ 2  3  0

julia> to_table(nplm)
6×3 Array{Any,2}:
 "a"  "a"  0
 "a"  "b"  1
 "a"  "c"  2
 "b"  "b"  0
 "b"  "c"  3
 "c"  "c"  0

```
