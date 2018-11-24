var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "PairwiseListMatrices",
    "title": "PairwiseListMatrices",
    "category": "page",
    "text": "CurrentModule = PairwiseListMatrices"
},

{
    "location": "#PairwiseListMatrices-1",
    "page": "PairwiseListMatrices",
    "title": "PairwiseListMatrices",
    "category": "section",
    "text": ""
},

{
    "location": "#Introduction-1",
    "page": "PairwiseListMatrices",
    "title": "Introduction",
    "category": "section",
    "text": "PairwiseListMatrices allows you to represent a (squared) symmetric matrix as a list of the values in the upper or lower triangular part of the matrix. Those matrices are common for representing pairwise measures/comparisons between the elements of one group when the used metric/distance satisfies the symmetry condition. Also the adjacency matrices of undirected graphs can be represented with this kind of list/matrices.  "
},

{
    "location": "#Installation-1",
    "page": "PairwiseListMatrices",
    "title": "Installation",
    "category": "section",
    "text": "This package is registered on Julia\'s METADATA.jl and it can be installed through the Julia\'s REPL:  Pkg.add(\"PairwiseListMatrices\")If the package is installed on your system, you can load it with:  using PairwiseListMatrices"
},

{
    "location": "#Simple-example-1",
    "page": "PairwiseListMatrices",
    "title": "Simple example",
    "category": "section",
    "text": "The following symmetric matrix has 9 values. Their values could be thought as pairwise measures between 3 elements:  matrix = [  0 10 20\n           10  0 30\n           20 30  0 ]Since all the diagonal members are zeros, this matrix could be represented as a vector/list of the 3 values on the triangular part:  list = [10, 20, 30]The type PairwiseListMatrix, defined in this module, can be used for working with the list as a full symmetric matrix.  using PairwiseListMatrices\nplm = PairwiseListMatrix(list)"
},

{
    "location": "#Implementation-1",
    "page": "PairwiseListMatrices",
    "title": "Implementation",
    "category": "section",
    "text": "If you are performing pairwise measures over N elements, storing all the N*N values of a Matrix{T} represents sizeof(T)*(N*N) bytes of memory. However, the lower and upper triangular parts of the matrix are identical and could be stored in a single list. In this way, you are storing the green value only once:  (Image: )The diagonal values should be stored, since they could change at any time (i.e. yellow value). So you need sizeof(T)*(N) bytes for storing the diagonal values on a vector and sizeof(T)*(N*(N-1))/2 bytes for storing the lower or upper triangular part of the matrix. The type PairwiseListMatrix{T, diagonal, VT} represents the symmetric matrix using only sizeof(T)*(N*(N+1))/2 bytes instead of sizeof(T)*(N*N) bytes, saving almost 50% of the memory (the percent depends on N):  using Plots\ngr()\nplot(      2:550,\n           N -> 100.0 - ( 100.0 * div(N*(N+1), 2) / (N*N) ),\n           xlab = \"N\",\n           ylab = \"% of saved memory\",\n           legend = nothing        )\npng(\"curve.png\") # hide\nnothing # hide(Image: )As you can see in the schematic diagram, the difference between PairwiseListMatrix{T, true, VT} and PairwiseListMatrix{T, false, VT} is where the diagonal values are stored. All PairwiseListMatrix{T, diagonal, VT} have a list field for storing the values. If diagonal is true, the diagonal values are included in the list (i.e. yellow value) and the diag vector is empty. But if the diagonal value is false the diagonal values are stored in the diag vector.  mutable struct PairwiseListMatrix{T,diagonal,VT} <: AbstractArray{T, 2}\n    list::VT\n    diag::VT\n    nelements::Int\n    ...\nendThe number of elements in the pairwise measure/comparisons or the number of nodes in the undirected graph is stored in nelements and used in indexing operations. This allows you to index the object like any other matrix.  The PairwiseListMatrix can be wrapped in a NamedArray (from the package NamedArrays) to allow the access of elements using labels. The function setlabel can be used to create this object easily. For example, using the matrix of the figure and storing the diagonal values in the list:using PairwiseListMatrices\nplm = PairwiseListMatrix([1,1,0,1,0,1,1,0,0,0], true)\nnplm = setlabels(plm, [\"A\",\"B\",\"C\",\"D\"])\nnplm[\"B\",\"C\"]You can also create the matrix with the list without the diagonal values and fill the diagonal values after that:  using PairwiseListMatrices\nplm = PairwiseListMatrix([1,0,1,1,1,0], false)\nnplm = setlabels(plm, [\"A\",\"B\",\"C\",\"D\"])\nnplm[\"A\",\"A\"] = 1\nnplm"
},

{
    "location": "#Ploting-1",
    "page": "PairwiseListMatrices",
    "title": "Ploting",
    "category": "section",
    "text": "You can use the Plots package to visualize this matrices quickly as heat maps. If you are looking for more complex visualization, you can use the PlotRecipes package. This last package provides arc diagram, chord diagram/circos and other graphplots (since those matrices could be a representation for an adjacency matrix/list of an undirected graph).  using PairwiseListMatrices\nplm = PairwiseListMatrix([1,1,0,1,0,1,1,0,0,0], true)\nnplm = setlabels(plm, [\"A\",\"B\",\"C\",\"D\"])\n\nusing Plots\ngr()\nplot(nplm)\npng(\"heatmap.png\") # hide\nnothing # hide(Image: )  "
},

{
    "location": "#Benchmark-1",
    "page": "PairwiseListMatrices",
    "title": "Benchmark",
    "category": "section",
    "text": "PairwiseListMatrix is faster than a full matrix to make operation like sum and mean in the whole matrix, since it is cache efficient. However it is slower than a full matrix for reducing along dimensions.Creation benchmark\nStatistics benchmark"
},

{
    "location": "api/#",
    "page": "API",
    "title": "API",
    "category": "page",
    "text": ""
},

{
    "location": "api/#API-1",
    "page": "API",
    "title": "API",
    "category": "section",
    "text": ""
},

{
    "location": "api/#PairwiseListMatrices.PairwiseListMatrix",
    "page": "API",
    "title": "PairwiseListMatrices.PairwiseListMatrix",
    "category": "type",
    "text": "PairwiseListMatrix{T, diagonal, VT} is a (squared) symmetric matrix that stores a list of type VT with values of type T for the pairwise comparison/evaluation of nelements. If diagonal is true the first element of the list is 1, 1 otherwise is 1, 2. If diagonal is false the diagonal values are stored in a vector on the diag field.\n\n\n\n\n\n"
},

{
    "location": "api/#Creation-1",
    "page": "API",
    "title": "Creation",
    "category": "section",
    "text": "PairwiseListMatrix"
},

{
    "location": "api/#PairwiseListMatrices.hasdiagonal",
    "page": "API",
    "title": "PairwiseListMatrices.hasdiagonal",
    "category": "function",
    "text": "Returns true if the list has diagonal values.\n\n\n\n\n\n"
},

{
    "location": "api/#PairwiseListMatrices.getlist",
    "page": "API",
    "title": "PairwiseListMatrices.getlist",
    "category": "function",
    "text": "Retuns the list vector.\n\n\n\n\n\n"
},

{
    "location": "api/#PairwiseListMatrices.getdiag",
    "page": "API",
    "title": "PairwiseListMatrices.getdiag",
    "category": "function",
    "text": "Retuns the diag vector (which contains the diagonal values if diagonal is false).\n\n\n\n\n\n"
},

{
    "location": "api/#Getters-1",
    "page": "API",
    "title": "Getters",
    "category": "section",
    "text": "hasdiagonal\ngetlist\ngetdiag"
},

{
    "location": "api/#PairwiseListMatrices.lengthlist",
    "page": "API",
    "title": "PairwiseListMatrices.lengthlist",
    "category": "function",
    "text": "Returns the length of the list field\n\n\n\n\n\nReturns the list length needed for a pairwise measures or comparisons of nelements. If diagonal is true, diagonal values are included in the list.\n\njulia> using PairwiseListMatrices\n\njulia> plm = PairwiseListMatrix([1, 2, 3, 4, 5, 6], false)\n4×4 PairwiseListMatrix{Int64,false,Array{Int64,1}}:\n 0  1  2  3\n 1  0  4  5\n 2  4  0  6\n 3  5  6  0\n\njulia> lengthlist(4, false)\n6\n\njulia> lengthlist(plm)\n6\n\n\n\n\n\n\n"
},

{
    "location": "api/#PairwiseListMatrices.ij2k",
    "page": "API",
    "title": "PairwiseListMatrices.ij2k",
    "category": "function",
    "text": "Returns the k index of the list from the indixes i and j with i<j from a matrix of nelements by nelements. diagonal should be true or Val{true} if the diagonal values are on the list. You must not use it with i>j.\n\njulia> using PairwiseListMatrices\n\njulia> plm = PairwiseListMatrix([10,20,30,40,50,60], true)\n3×3 PairwiseListMatrix{Int64,true,Array{Int64,1}}:\n 10  20  30\n 20  40  50\n 30  50  60\n\njulia> ij2k(1, 2, 3, true)\n2\n\njulia> getlist(plm)[2]\n20\n\n\n\n\n\n\n"
},

{
    "location": "api/#PairwiseListMatrices.diagonal",
    "page": "API",
    "title": "PairwiseListMatrices.diagonal",
    "category": "function",
    "text": "Returns a vector of type VT from a PairwiseListMatrix{T, false, VT} that has the diagonal values.\n\njulia> using PairwiseListMatrices\n\njulia> plm = PairwiseListMatrix([10,20,30,40,50,60], true)\n3×3 PairwiseListMatrix{Int64,true,Array{Int64,1}}:\n 10  20  30\n 20  40  50\n 30  50  60\n\njulia> diagonal(plm)\n3-element Array{Int64,1}:\n 10\n 40\n 60\n\n\n\n\n\n\n"
},

{
    "location": "api/#Helpers-1",
    "page": "API",
    "title": "Helpers",
    "category": "section",
    "text": "lengthlist\nij2k\ndiagonal"
},

{
    "location": "api/#PairwiseListMatrices.@iteratelist",
    "page": "API",
    "title": "PairwiseListMatrices.@iteratelist",
    "category": "macro",
    "text": "The macro @iteratelist writes a for loop over the list but avoiding getfield calls inside the loop. The first argument of the macro is the PairwiseListMatrix that is going to be iterated and the second is the body of the loop. In the body list will be the list field of the PairwiseListMatrix and k the index over that list. Other variables should be interpolated in a quote. You must not modify the value of k.\n\njulia> using PairwiseListMatrices\n\njulia> PLM = PairwiseListMatrix([1,2,3], false)\n3×3 PairwiseListMatrix{Int64,false,Array{Int64,1}}:\n 0  1  2\n 1  0  3\n 2  3  0\n\njulia> @iteratelist PLM println(list[k])\n1\n2\n3\n\n\n\n\n\n\n"
},

{
    "location": "api/#PairwiseListMatrices.@iteratediag",
    "page": "API",
    "title": "PairwiseListMatrices.@iteratediag",
    "category": "macro",
    "text": "The macro @iteratediag writes a for loop over the diag field of a PairwiseListMatrix{T,false,VT} but avoiding calls to getfield inside the loop. The first argument of the macro is the PairwiseListMatrix that is going to be iterated and the second is the body of the loop. In the body diag will be the diag field of the PairwiseListMatrix and k the index over that vector. Other variables should be interpolated in a quote. You must not modify the value of k.\n\njulia> using PairwiseListMatrices\n\njulia> PLM = PairwiseListMatrix([1,2,3], false)\n3×3 PairwiseListMatrix{Int64,false,Array{Int64,1}}:\n 0  1  2\n 1  0  3\n 2  3  0\n\njulia> @iteratediag PLM diag[k] += 10k\n\njulia> PLM\n3×3 PairwiseListMatrix{Int64,false,Array{Int64,1}}:\n 10   1   2\n  1  20   3\n  2   3  30\n\n\n\n\n\n\n"
},

{
    "location": "api/#PairwiseListMatrices.@iterateupper",
    "page": "API",
    "title": "PairwiseListMatrices.@iterateupper",
    "category": "macro",
    "text": "The macro @iterateupper iterates over the upper triangular part of the PairwiseListMatrix that is given as first argument. The second argument should be true if the diagonal need to be included in the iteration or false otherwise. The last argument is the body of the loop, where list is the list and diag fields of the PairwiseListMatrix and k is the index over that list. You can also use the respective i and j indexes for that position k in the upper triangular part of the matrix. Other variables should be interpolated in a quote. You must not modify the values of i, j or k.\n\njulia> using PairwiseListMatrices\n\njulia> PLM = PairwiseListMatrix([1,2,3], true)\n2×2 PairwiseListMatrix{Int64,true,Array{Int64,1}}:\n 1  2\n 2  3\n\njulia> mat = zeros(Int, 2, 2)\n2×2 Array{Int64,2}:\n 0  0\n 0  0\n\njulia> let mat = mat # To avoid using global\n           @iterateupper PLM true :($mat)[i,j] = list[k]\n       end\n\njulia> mat\n2×2 Array{Int64,2}:\n 1  2\n 0  3\n\n\n\n\n\n\n"
},

{
    "location": "api/#Macros-1",
    "page": "API",
    "title": "Macros",
    "category": "section",
    "text": "@iteratelist\n@iteratediag\n@iterateupper"
},

{
    "location": "api/#PairwiseListMatrices.getlabels",
    "page": "API",
    "title": "PairwiseListMatrices.getlabels",
    "category": "function",
    "text": "It gets the labels of a PairwiseListMatrix.\n\njulia> using PairwiseListMatrices\n\njulia> plm  = PairwiseListMatrix(ones(3), false)\n3×3 PairwiseListMatrix{Float64,false,Array{Float64,1}}:\n 0.0  1.0  1.0\n 1.0  0.0  1.0\n 1.0  1.0  0.0\n\njulia> getlabels(plm)\n3-element Array{String,1}:\n \"1\"\n \"2\"\n \"3\"\n\njulia> nplm  = setlabels(plm, [\"a\",\"b\",\"c\"])\n3×3 Named PairwiseListMatrix{Float64,false,Array{Float64,1}}\nA ╲ B │   a    b    c\n──────┼──────────────\na     │ 0.0  1.0  1.0\nb     │ 1.0  0.0  1.0\nc     │ 1.0  1.0  0.0\n\njulia> getlabels(nplm)\n3-element Array{String,1}:\n \"a\"\n \"b\"\n \"c\"\n\n\n\n\n\n\n"
},

{
    "location": "api/#PairwiseListMatrices.setlabels",
    "page": "API",
    "title": "PairwiseListMatrices.setlabels",
    "category": "function",
    "text": "Creates a Named PairwiseListMatrix.\n\njulia> using PairwiseListMatrices\n\njulia> plm  = PairwiseListMatrix(ones(3), false)\n3×3 PairwiseListMatrix{Float64,false,Array{Float64,1}}:\n 0.0  1.0  1.0\n 1.0  0.0  1.0\n 1.0  1.0  0.0\n\njulia> nplm  = setlabels(plm, [\"a\",\"b\",\"c\"])\n3×3 Named PairwiseListMatrix{Float64,false,Array{Float64,1}}\nA ╲ B │   a    b    c\n──────┼──────────────\na     │ 0.0  1.0  1.0\nb     │ 1.0  0.0  1.0\nc     │ 1.0  1.0  0.0\n\n\n\n\n\n\n"
},

{
    "location": "api/#PairwiseListMatrices.setlabels!",
    "page": "API",
    "title": "PairwiseListMatrices.setlabels!",
    "category": "function",
    "text": "It changes the labels of a Named PairwiseListMatrix\n\n\n\n\n\n"
},

{
    "location": "api/#Labels/Names-1",
    "page": "API",
    "title": "Labels/Names",
    "category": "section",
    "text": "getlabels\nsetlabels\nsetlabels!"
},

{
    "location": "api/#Base.join",
    "page": "API",
    "title": "Base.join",
    "category": "function",
    "text": "This function join two PairwiseListMatrices by their labels, returning two PairwiseListMatrices with same size and labels. There are 4 kinds of joins:\n\n:inner : Intersect. The output matrices only include the labels that are in both PairwiseListMatrices\n:outer : Union. Include the labels of the two PairwiseListMatrices.\n:left : Only use labels from the first argument.\n:right : Only use labels from the second argument.\n\nNaNs are filled in where needed to complete joins. The default value for missing values can be changed passing a tuple to missing.\n\njulia> using PairwiseListMatrices\n\njulia> l = setlabels(PairwiseListMatrix([1.,2.,3.], false), [\"a\",\"b\",\"c\"]) # a b c\n3×3 Named PairwiseListMatrix{Float64,false,Array{Float64,1}}\nA ╲ B │   a    b    c\n──────┼──────────────\na     │ 0.0  1.0  2.0\nb     │ 1.0  0.0  3.0\nc     │ 2.0  3.0  0.0\n\njulia> r = setlabels(PairwiseListMatrix([1.,2.,3.], false), [\"b\",\"c\",\"d\"]) # b c d\n3×3 Named PairwiseListMatrix{Float64,false,Array{Float64,1}}\nA ╲ B │   b    c    d\n──────┼──────────────\nb     │ 0.0  1.0  2.0\nc     │ 1.0  0.0  3.0\nd     │ 2.0  3.0  0.0\n\njulia> join(l, r, kind=:inner) # b c\n(2×2 Named PairwiseListMatrix{Float64,false,Array{Float64,1}}\nA ╲ B │   b    c\n──────┼─────────\nb     │ 0.0  3.0\nc     │ 3.0  0.0, 2×2 Named PairwiseListMatrix{Float64,false,Array{Float64,1}}\nA ╲ B │   b    c\n──────┼─────────\nb     │ 0.0  1.0\nc     │ 1.0  0.0)\n\njulia> join(l, r, kind=:outer) # a b c d\n(4×4 Named PairwiseListMatrix{Float64,false,Array{Float64,1}}\nA ╲ B │   a    b    c    d\n──────┼───────────────────\na     │ 0.0  1.0  2.0  NaN\nb     │ 1.0  0.0  3.0  NaN\nc     │ 2.0  3.0  0.0  NaN\nd     │ NaN  NaN  NaN  NaN, 4×4 Named PairwiseListMatrix{Float64,false,Array{Float64,1}}\nA ╲ B │   a    b    c    d\n──────┼───────────────────\na     │ NaN  NaN  NaN  NaN\nb     │ NaN  0.0  1.0  2.0\nc     │ NaN  1.0  0.0  3.0\nd     │ NaN  2.0  3.0  0.0)\n\n\n\n\n\n\n"
},

{
    "location": "api/#Join-1",
    "page": "API",
    "title": "Join",
    "category": "section",
    "text": "join"
},

{
    "location": "api/#PairwiseListMatrices.sum_nodiag",
    "page": "API",
    "title": "PairwiseListMatrices.sum_nodiag",
    "category": "function",
    "text": "Sum the values outside the diagonal\n\n\n\n\n\n"
},

{
    "location": "api/#PairwiseListMatrices.mean_nodiag",
    "page": "API",
    "title": "PairwiseListMatrices.mean_nodiag",
    "category": "function",
    "text": "Mean of the values outside the diagonal\n\n\n\n\n\n"
},

{
    "location": "api/#PairwiseListMatrices.zscore",
    "page": "API",
    "title": "PairwiseListMatrices.zscore",
    "category": "function",
    "text": "It\'s like zscore! but without modifying the PairwiseListMatrix\n\n\n\n\n\n"
},

{
    "location": "api/#PairwiseListMatrices.zscore!",
    "page": "API",
    "title": "PairwiseListMatrices.zscore!",
    "category": "function",
    "text": "This function takes a vector of PairwiseListMatrix objects and a PairwiseListMatrix and fill the matrix with the zscore value using the median and std of the vector.\n\n\n\n\n\n"
},

{
    "location": "api/#Statistics-1",
    "page": "API",
    "title": "Statistics",
    "category": "section",
    "text": "sum_nodiag\nmean_nodiag\nzscore\nzscore!"
},

{
    "location": "api/#PairwiseListMatrices.to_table",
    "page": "API",
    "title": "PairwiseListMatrices.to_table",
    "category": "function",
    "text": "Creates a Matrix{Any}, labels are stored in the columns 1 and 2, and the values in the column 3. Diagonal values are included by default.\n\njulia> using PairwiseListMatrices\n\njulia> plm = PairwiseListMatrix([10,20,30], false)\n3×3 PairwiseListMatrix{Int64,false,Array{Int64,1}}:\n  0  10  20\n 10   0  30\n 20  30   0\n\njulia> to_table(plm)\n6×3 Array{Any,2}:\n \"1\"  \"1\"   0\n \"1\"  \"2\"  10\n \"1\"  \"3\"  20\n \"2\"  \"2\"   0\n \"2\"  \"3\"  30\n \"3\"  \"3\"   0\n\njulia> to_table(plm, diagonal=false)\n3×3 Array{Any,2}:\n \"1\"  \"2\"  10\n \"1\"  \"3\"  20\n \"2\"  \"3\"  30\n\n\n\n\n\n\n"
},

{
    "location": "api/#PairwiseListMatrices.to_dict",
    "page": "API",
    "title": "PairwiseListMatrices.to_dict",
    "category": "function",
    "text": "It takes a PairwiseListMatrix and converts it to a Dict of Symbols to arrays. The returned dictionary can be easily converted into a DataFrame.\n\njulia> using PairwiseListMatrices, DataFrames\n\njulia> nplm = setlabels(PairwiseListMatrix([10,20,30], false), [\"a\",\"b\",\"c\"])\n3×3 Named PairwiseListMatrix{Int64,false,Array{Int64,1}}\nA ╲ B │  a   b   c\n──────┼───────────\na     │  0  10  20\nb     │ 10   0  30\nc     │ 20  30   0\n\njulia> dict = to_dict(nplm, diagonal=false)\nDict{Symbol,Array{T,1} where T} with 3 entries:\n  :values => [10, 20, 30]\n  :j      => [\"b\", \"c\", \"c\"]\n  :i      => [\"a\", \"a\", \"b\"]\n\njulia> DataFrame(dict)\n3×3 DataFrames.DataFrame\n│ Row │ i      │ j      │ values │\n│     │ String │ String │ Int64  │\n├─────┼────────┼────────┼────────┤\n│ 1   │ a      │ b      │ 10     │\n│ 2   │ a      │ c      │ 20     │\n│ 3   │ b      │ c      │ 30     │\n\n\n\n\n\n"
},

{
    "location": "api/#PairwiseListMatrices.from_table",
    "page": "API",
    "title": "PairwiseListMatrices.from_table",
    "category": "function",
    "text": "Creation of a PairwiseListMatrix from a Matrix, DataFrame or similar structure. By default the columns with the labels for i (slow) and j (fast) are 1 and 2. Values are taken from the column 3 by default.\n\njulia> using PairwiseListMatrices, Pkg, DelimitedFiles\n\njulia> import PairwiseListMatrices\n\njulia> filename = joinpath(dirname(pathof(PairwiseListMatrices)), \"..\", \"test\", \"example.csv\");\n\njulia> dat = readdlm(filename, \',\')\n3×3 Array{Any,2}:\n \"A\"  \"B\"  10\n \"A\"  \"C\"  20\n \"B\"  \"C\"  30\n\njulia> from_table(dat, false)\n3×3 Named PairwiseListMatrix{Any,false,Array{Any,1}}\nA ╲ B │       A        B        C\n──────┼──────────────────────────\nA     │ nothing       10       20\nB     │      10  nothing       30\nC     │      20       30  nothing\n\nThis is also useful to create a PairwiseListMatrix from a DataFrame:\n\njulia> using PairwiseListMatrices, DataFrames, CSV, Pkg\n\njulia> import PairwiseListMatrices\n\njulia> filename = joinpath(dirname(pathof(PairwiseListMatrices)), \"..\", \"test\", \"example.csv\");\n\njulia> df = CSV.read(filename, header=false)\n3×3 DataFrames.DataFrame\n│ Row │ Column1 │ Column2 │ Column3 │\n│     │ String⍰ │ String⍰ │ Int64⍰  │\n├─────┼─────────┼─────────┼─────────┤\n│ 1   │ A       │ B       │ 10      │\n│ 2   │ A       │ C       │ 20      │\n│ 3   │ B       │ C       │ 30      │\n\njulia> from_table(df, false)\n3×3 Named PairwiseListMatrix{Union{Missing, Int64},false,Array{Union{Missing, Int64},1}}\nA ╲ B │       A        B        C\n──────┼──────────────────────────\nA     │ missing       10       20\nB     │      10  missing       30\nC     │      20       30  missing\n\n\n\n\n\n"
},

{
    "location": "api/#DelimitedFiles.writedlm",
    "page": "API",
    "title": "DelimitedFiles.writedlm",
    "category": "function",
    "text": "This function takes the filename as first argument and a PairwiseListMatrix as second argument. If the diagonal keyword argument is true (default), the diagonal is included in the output. The keyword argument delim (by default is \'	\') allows to modified the character used as delimiter.\n\njulia> using PairwiseListMatrices, DelimitedFiles\n\njulia> plm  = PairwiseListMatrix(trues(3), false)\n3×3 PairwiseListMatrix{Bool,false,BitArray{1}}:\n false   true   true\n  true  false   true\n  true   true  false\n\njulia> writedlm(\"example.csv\", plm, diagonal=false, delim=\',\')\n\njulia> println(read(\"example.csv\", String))\n1,2,true\n1,3,true\n2,3,true\n\n\n\n\n\n\n"
},

{
    "location": "api/#IO-1",
    "page": "API",
    "title": "IO",
    "category": "section",
    "text": "to_table\nto_dict\nfrom_table\nwritedlm"
},

]}
