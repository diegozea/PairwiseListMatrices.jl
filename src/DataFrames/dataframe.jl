export to_dataframe, from_dataframe

"""
Creates a `DataFrame` from a `PairwiseListMatrix`. The second positional argument `diag` is `true` by default
and the diagonal values are included in the `DataFrame`.
Labels are stored in the columns 1 (`i`) and 2 (`j`), and the values in the column 3 with `colname` as `:value` by default.

```
julia> using DataFrames, PairwiseListMatrices

julia> PLM  = PairwiseListMatrix([1,2,3,4,5,6], true)
3x3 PairwiseListMatrices.PairwiseListMatrix{Int64,true}:
 1  2  3
 2  4  5
 3  5  6

julia> to_dataframe(PLM, false, colname=:score)
3x3 DataFrames.DataFrame
| Row | i | j | score |
|-----|---|---|-------|
| 1   | 1 | 2 | 2     |
| 2   | 1 | 3 | 3     |
| 3   | 2 | 3 | 5     |

```
"""
function to_dataframe{T,D,TV}(plm::PairwiseListMatrix{T,D,TV},
                              diag::Bool=true;
                              labels::Vector{String} = labels(plm),
                              colname::Symbol=:value)
    N = plm.nelements
    df = DataFrames.DataFrame([String, String, T],
                              [:i, :j, colname],
                              diag ? div(N*(N+1),2) : div(N*(N-1),2))
    t = 0
    @iterateupper plm diag begin
        t += 1
        df[t, :i] = labels[i]
        df[t, :j] = labels[j]
        df[t, colname] = list[k]
    end
    df
end

"""
Creation of a `PairwiseListMatrix` from a `DataFrame`.
By default the columns with the labels for i (slow) and j (fast) are 1 and 2. Values are taken from the column 3 by default.
The second argument `diagonal` should be true if the diagonal values are included in the `DataFrame`.
"""
function from_dataframe(df::DataFrames.DataFrame,
                        diagonal::Bool;
                        labelcols::Vector{Int} = [1,2],
                        valuecol::Int = 3)
    plm = PairwiseListMatrix(df[:, valuecol], diagonal)
    nplm = NamedArray(plm)
    if length(labelcols) == 2
        labels = unique(vcat(df[labelcols[1]], df[labelcols[2]]))
        setlabels!(nplm, String[ string(lab) for lab in labels ])
    end
    nplm
end
