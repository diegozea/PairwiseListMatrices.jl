## PairwiseListMatrices.jl Release Notes

## Changes from v0.6 to v0.7

PairwiseListMatrices v0.7 requires Julia 0.7/1.0.
PairwiseListMatrices v0.6 was the last version with Julia 0.6 support.

* `writecsv` changed to `writedlm` using `delim=','`.   

* A type parameter to indicate the wrapped matrix type was added to
`NamedResidueMatrix`, which is now `NamedResidueMatrix{Array{Residue,2}}`.  

### Changes from v0.5 to v0.6

PairwiseListMatrices v0.6 requires Julia v0.6.
PairwiseListMatrices v0.5 was the last version with Julia 0.5 support.

### Changes from v0.4 to v0.5

* It has better support for Julia 0.5, e.g. faster `map` and `broadcast`.  

* The `label` field of `PairwiseListMatrix` is eliminated, `NamedArrays` are used instead
for storing element names.  

* `PairwiseListMatrix` now takes the type of the list vector as type parameter. This gives
better support for different vector types e.g. `BitArray`.   

* `RecipesBase` are used for `plot` the matrices `using Plots`.  

* `from_` and `to_dataframe` are deprecated in favor of a more general `from_table` and
`DataFrame(to_dict(...))` respectively.  

* It solves some bugs.  

### Changes from v0.3 to v0.4

PairwiseListMatrices v0.4 requires Julia v0.5.
PairwiseListMatrices v0.3 was the last version with Julia 0.4 support.

### Changes from v0.2 to v0.3

* `join` two PairwiseListMatrices by their labels

* `PairwiseListMatrix` nows works with any `AbstractVector`. i.e. `DataArray`s

### Changes from v0.1 to v0.2

* `hasdiagonal` returns the `diagonal` boolean value of `PairwiseListMatrix{T, diagonal}`

* `getlist` and `getdiag` are defined for getting the `list`or `diag` fields of a `PairwiseListMatrix`.

* `lengthlist` was added for returning the length of a list given a `PairwiseListMatrix` or a number of elements.

* `ij2k` was added for getting the indexes on the list given the `i` and `j` on the upper triangular part of the matrix.

* `@iterateupper`, `@iteratelist` and `@iteratediag` for easy iterations over a `PairwiseListMatrix` and its fields.

* `protovis` was added for plotting a `PairwiseListMatrix` as a network using the javascript protovis library.

* `writecsv` and `writedlm` for `PairwiseListMatrix`

* `from_table` and `to_table` go back and forth between `PairwiseListMatrix` and a three columns `Matrix`.

* `from_dataframe` and `to_dataframe` go back and forth between `PairwiseListMatrix` and a three columns `DataFrame`.

* `zscore` and `zscore!` were added and works on `Vector{PairwiseListMatrix{T, diagonal}}`

* `triu` is now defined for `PairwiseListMatrix`

* `convert` goes back and forth between `Matrix{T}` and `PairwiseListMatrix`

* `./` is now defined for PairwiseListMatrix of integers and the results are always `Float64`.

* The type parameter `L`, which had been indicating the type of the labels, was eliminated:  `PairwiseListMatrix{T, L, diagonal}` is now `PairwiseListMatrix{T, diagonal}`.

* Package precompilation is now allowed

PairwiseListMatrices v0.2 also includes several bug fixes.
