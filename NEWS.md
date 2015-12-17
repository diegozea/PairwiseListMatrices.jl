## PairwiseListMatrices.jl Release Notes

### Changes from v0.1 to v0.2

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
