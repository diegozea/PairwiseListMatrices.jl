using PairwiseListMatrices
using DataFrames
using NamedArrays
using Base.Test

@testset "PairwiseListMatrices" begin

    list = [1,-2,3]
    list_diag_triu =  [ 1 -2
                        0 3 ]
    list_diag_tril =  [ 1 0
                       -2 3 ]
    list_triu = [ 0 1 -2
                  0 0 3
                  0 0 0 ]
    list_tril = [ 0 0 0
                  1 0 0
                 -2 3 0 ]
    list_diag_sym = [ 1 -2
                     -2 3 ]
    list_sym = [ 0 1 -2
                 1 0 3
                -2 3 0 ]

    pld_triu = UpperTriangular(PairwiseListMatrix(list, true))
    pld_tril = LowerTriangular(PairwiseListMatrix(list, true))
    pl_triu  = UpperTriangular(PairwiseListMatrix(list))
    pl_tril  = LowerTriangular(PairwiseListMatrix(list))
    pld_sym  = PairwiseListMatrix(list, true)
    pl_sym   = PairwiseListMatrix(list)

    @testset "Creation from list" begin

        @test pld_triu == list_diag_triu
        @test pld_tril == list_diag_tril
        @test pl_triu  == list_triu
        @test pl_tril  == list_tril
        @test pld_sym  == list_diag_sym
        @test pl_sym   == list_sym
    end

    @testset "Getters" begin

        @test getlist(pld_sym) == pld_sym.list
        @test getlist(pl_sym)  ==  pl_sym.list

        @test getdiag(pld_sym) == pld_sym.diag
        @test getdiag(pl_sym)  ==  pl_sym.diag
    end

    @testset "Length" begin

        @test lengthlist(2, true)  == lengthlist(pld_sym)
        @test lengthlist(3, false) == lengthlist(pl_sym)
    end

    @testset "Convert" begin

        @test list_diag_sym ==
                convert(PairwiseListMatrix{Int, true, Vector{Int}},  list_diag_sym)
        @test list_diag_sym ==
                convert(PairwiseListMatrix{Int, false, Vector{Int}}, list_diag_sym)

        @test convert(Matrix{Float64}, list_sym) ==
                convert(PairwiseListMatrix{Float64, true, Vector{Float64}},  list_sym)
        @test convert(Matrix{Float64}, list_sym) ==
                convert(PairwiseListMatrix{Float64, false, Vector{Float64}}, list_sym)

        @test convert(Matrix{Int}, pl_sym)  == list_sym
        @test convert(Matrix{Int}, pld_sym) == list_diag_sym

        @test convert(PairwiseListMatrix{Int, false, Vector{Int}}, pld_sym) == pld_sym
        @test convert(PairwiseListMatrix{Int, true, Vector{Int}}, pl_sym)   == pl_sym
    end

    @testset "Unary operations" begin

        for f in [ -, abs, x -> sqrt(abs(x))]
            @test f(pld_triu) == f(list_diag_triu)
            @test f(pld_tril) == f(list_diag_tril)
            @test f(pl_triu ) == f(list_triu)
            @test f(pl_tril ) == f(list_tril)
            @test f(pld_sym ) == f(list_diag_sym)
            @test f(pl_sym  ) == f(list_sym)
        end
    end

    @testset "Binary operations" begin

        for f in [ -, +, .-, .+, .* ]
            @test f(pld_triu, pld_triu) == f(list_diag_triu, list_diag_triu)
            @test f(pld_tril, pld_tril) == f(list_diag_tril, list_diag_tril)
            @test f(pl_triu , pl_triu ) == f(list_triu,      list_triu     )
            @test f(pl_tril , pl_tril ) == f(list_tril,      list_tril     )
            @test f(pld_sym , pld_sym ) == f(list_diag_sym,  list_diag_sym )
            @test f(pl_sym  , pl_sym  ) == f(list_sym,       list_sym      )
        end

        @test pld_sym ./ pld_sym == list_diag_sym ./ list_diag_sym

        # https://github.com/JuliaLang/julia/issues/19615
        for f in [ .*, ./, / ]
            @test f(pld_triu, 2) == f(list_diag_triu, 2)
            @test f(pld_tril, 2) == f(list_diag_tril, 2)
            @test f(pl_triu , 2) == f(list_triu,      2)
            @test f(pl_tril , 2) == f(list_tril,      2)
        end

        for f in [ -, +, .-, .+, .*, ./, / ]
            @test f(pld_sym , 2) == f(list_diag_sym,  2)
            @test f(pl_sym  , 2) == f(list_sym,       2)
        end

        @testset "./ with Float64 & Int" begin

            for value in  (4, 4.0)
                list = PairwiseListMatrix([.5, .4, .3])
                result =  list ./ value
                @test isa(result, PairwiseListMatrix{Float64,false,Vector{Float64}})
                @test sum(result) ≈ 0.6

                list = PairwiseListMatrix([.5, .4, .3], true)
                result =  list ./ value
                @test isa(result, PairwiseListMatrix{Float64,true,Vector{Float64}})
                @test sum(result) ≈ (0.5 + 0.4*2.0 + 0.3)/4
            end
        end
    end

    @testset "Transpose" begin

        @test transpose(pld_triu) == pld_tril == list_diag_triu.'
        @test transpose(pld_tril) == pld_triu == list_diag_tril.'
        @test transpose(pl_triu) == pl_tril == list_triu.'
        @test transpose(pl_tril) == pl_triu == list_tril.'

        @test ctranspose(pld_triu) == pld_tril == list_diag_triu'
        @test ctranspose(pld_tril) == pld_triu == list_diag_tril'
        @test ctranspose(pl_triu) == pl_tril == list_triu'
        @test ctranspose(pl_tril) == pl_triu == list_tril'

        @test transpose(pld_sym) == pld_sym == list_diag_sym.'
        @test transpose(pl_sym) == pl_sym == list_sym.'

        @test ctranspose(pld_sym) == pld_sym == list_diag_sym'
        @test ctranspose(pl_sym) == pl_sym == list_sym'

        @test transpose!(pl_sym) == pl_sym
        @test ctranspose!(pl_sym) == pl_sym
    end

    @testset "Linear algebra" begin

        @test pld_triu * pld_triu == list_diag_triu * list_diag_triu
        @test pld_tril * pld_tril == list_diag_tril * list_diag_tril
        @test pl_triu  * pl_triu  == list_triu      * list_triu
        @test pl_tril  * pl_tril  == list_tril      * list_tril
        @test pld_sym  * pld_sym  == list_diag_sym  * list_diag_sym
        @test pl_sym   * pl_sym   == list_sym       * list_sym
        @test Symmetric(pl_sym) * Symmetric(pl_sym) == Symmetric(list_sym) * Symmetric(list_sym)

        @test pld_triu / pld_triu == list_diag_triu / list_diag_triu
        @test pld_tril / pld_tril == list_diag_tril / list_diag_tril
        @test pld_sym  / pld_sym  == list_diag_sym  / list_diag_sym
        @test pl_sym   / pl_sym   == list_sym       / list_sym

        @test svd(pld_triu) == svd(list_diag_triu)
        @test svd(pld_tril) == svd(list_diag_tril)
        @test svd(pl_triu ) == svd(list_triu)
        @test svd(pl_tril ) == svd(list_tril)
        @test svd(pld_sym ) == svd(list_diag_sym)
        @test svd(pl_sym  ) == svd(list_sym)
    end

    @testset "Stats" begin

        for f in [ mean, std ]
            @test f(pld_triu) == f(list_diag_triu)
            @test f(pld_tril) == f(list_diag_tril)
            @test f(pl_triu ) == f(list_triu)
            @test f(pl_tril ) == f(list_tril)
            @test f(pld_sym ) == f(list_diag_sym)
            @test f(pl_sym  ) == f(list_sym)

            @test f(pld_sym,1) == f(list_diag_sym, 1)
            @test f(pl_sym ,1) == f(list_sym, 1)
            @test f(pld_sym,2) == f(list_diag_sym, 2)
            @test f(pl_sym ,2) == f(list_sym, 2)
        end

        @test cor(pld_sym ) == cor(list_diag_sym)
        @test cor(pl_sym  ) == cor(list_sym)
    end

    @testset "triu" begin

        @test triu(pld_sym ) == triu(list_diag_sym)
        @test triu(pl_sym  ) == triu(list_sym)

        @test triu(pld_sym,1) == triu(list_diag_sym, 1)
        @test triu(pl_sym, 1) == triu(list_sym, 1)

        @test_throws ErrorException triu!(pl_sym)
        @test_throws ErrorException triu!(pl_sym, 1)
    end

    # + and - definitions are ambiguous with:
    # +(A::AbstractArray{Bool,N<:Any}, x::Bool)
    #
    # @testset "Boolean binary op" begin
    #
    #     list = PairwiseListMatrix([true, true, true])
    #
    #     @test list - true == [ -1  0  0
    #                                  0 -1  0
    #                                  0  0 -1 ]
    #
    #     @test list + true == [ 1  2  2
    #                                 2  1  2
    #                                 2  2  1 ]
    # end

    @testset "Mean without diagonal" begin

        diag_false = PairwiseListMatrix([10,20,30], false)
        diag_true  = PairwiseListMatrix([0,10,20,0,30,0], true)

        for list in [diag_false, diag_true]
            @test vec( mean_nodiag(list, 1) ) == [15., 20., 25.]
            @test vec( mean_nodiag(list, 2) ) == [15., 20., 25.]
            @test mean_nodiag(list) == 20.
        end
    end
end

@testset "Named PairwiseListMatrix" begin

    list =  [1,2,3]
    labels= ["a", "b", "c"]
    labels_diag = ["A", "B"]

    list_diag = setlabels(PairwiseListMatrix(list, true), labels_diag)
    list = setlabels(PairwiseListMatrix(list), labels)

    @testset "set and get labels" begin

        @test isa(list_diag, NamedArray)
        @test isa(list, NamedArray)
        @test isa(list_diag.array, PairwiseListMatrix)
        @test isa(list.array, PairwiseListMatrix)

        @test_throws AssertionError setlabels!(list_diag, labels)
        @test_throws AssertionError setlabels!(list, labels_diag)

        @test getlabels(list_diag) == labels_diag
        @test getlabels(list) == labels

        setlabels!(list, ["A","B","C"])
        @test getlabels(list) == ["A","B","C"]
        setlabels!(list, ["a","b","c"])
    end

    @testset "get/setindex" begin

        for i in 1:2
            for j in 2:3
                @test list[labels[i], labels[j]] == list[labels[j], labels[i]]
            end
        end

        @test list_diag["A", "B"] == list_diag["B", "A"]

        list_diag["A", "B"] = 20
        @test list_diag["B", "A"] == 20
        list_diag["B", "A"] = 2
        @test list_diag["A", "B"] == 2
    end

    @testset "join" begin

        left = setlabels(PairwiseListMatrix(collect(1.:15.), true),
            String["a","b","c","d","e"])
        right = setlabels(PairwiseListMatrix(collect(1.:10.), false),
            String["a","e","i","o","u"])

        a, b = join(left, right)
        @test (a, b) == join(left, right, kind=:inner)
        @test getlabels(a) == getlabels(b)
        @test getlabels(a) == ["a", "e"]
        @test size(a) == (2,2)
        @test size(b) == (2,2)

        a, b = join(left, right, kind=:left)
        @test getlabels(a) == getlabels(b)
        @test getlabels(a) == getlabels(left)
        @test a == left
        @test b != right

        a, b = join(left, right, kind=:right)
        @test getlabels(a) == getlabels(b)
        @test getlabels(b) == getlabels(right)
        @test a != left
        @test b == right

        a, b = join(left, right, kind=:outer)
        @test getlabels(a) == getlabels(b)
        @test getlabels(a) == ["a", "b", "c", "d", "e", "i", "o", "u"]
        @test size(a) == (8,8)
        @test size(b) == (8,8)
        @test a != left
        @test b != right
    end
end

@testset "Vector{PairwiseListMatrix}" begin

    list_samples = [ PairwiseListMatrix(rand(10), true) for i in 1:100 ]
    list_diag_samples = [ PairwiseListMatrix(rand(6), false) for i in 1:100 ]
    full_samples = Matrix{Float64}[ full(mat) for mat in list_samples ]
    full_diag_samples = Matrix{Float64}[ full(mat) for mat in list_diag_samples ]

    for f in [sum, mean]
        @test f(list_samples) == f(full_samples)
        @test f(list_diag_samples) == f(full_diag_samples)
    end

    @test sum(list_samples[1:1]) == list_samples[1]
    @test sum(list_diag_samples[1:1]) == list_diag_samples[1]

    @test_throws ErrorException sum(PairwiseListMatrix{Int, true, Array{Int,1}}[])
    @test_throws ErrorException sum(PairwiseListMatrix{Int, true, Array{Int,1}}[
                                        PairwiseListMatrix([1, 1, 1], true),
                                        PairwiseListMatrix([1, 1, 1, 1, 1, 1], true) ])
    @test_throws ErrorException std(PairwiseListMatrix{Int, true, Array{Int,1}}[
                                        PairwiseListMatrix([1, 1, 1], true) ])

    full_samples = reshape(hcat(full_samples...), (4,4,100))
    full_diag_samples = reshape(hcat(full_diag_samples...), (4,4,100))
    # @test std(list_samples) ≈ std(full_samples,3)[:,:,1] # It uses n - 1
    @test std(list_samples) ≈ sqrt(mean(abs2.(full_samples .- mean(full_samples,3)),3))
    # @test std(list_diag_samples) ≈ std(full_diag_samples,3)[:,:,1] # It uses n - 1
    @test std(list_diag_samples) ≈ sqrt(mean(abs2.(full_diag_samples .- mean(full_diag_samples,3)),3))

    @testset "Z score" begin

        list = PairwiseListMatrix{Float64, false, Vector{Float64}}[
            PairwiseListMatrix([1., 1., 1.]),
            PairwiseListMatrix([3., 3., 3.]) ]
        mat = PairwiseListMatrix([2., 2., 2.])

        @test std(list) == PairwiseListMatrix([1., 1., 1.])
        @test mean(list) == mat
        @test PairwiseListMatrices.zscore(list, mat) == PairwiseListMatrix([0., 0., 0.])

        @test_throws MethodError PairwiseListMatrices.zscore(list, PairwiseListMatrix([2.,2.,2.], true))
        @test_throws ErrorException PairwiseListMatrices.zscore(list, PairwiseListMatrix([2.,2.,2.,2.,2.,2.]))

        list = PairwiseListMatrix{Float64, true, Vector{Float64}}[
            PairwiseListMatrix([1., 1., 1.], true),
            PairwiseListMatrix([3., 3., 3.], true) ]
        mat = PairwiseListMatrix([2., 2., 2.], true)

        @test std(list) == PairwiseListMatrix([1., 1., 1.], true)
        @test mean(list) == mat
        @test PairwiseListMatrices.zscore(list, mat) == PairwiseListMatrix([0., 0., 0.], true)

        @test_throws MethodError PairwiseListMatrices.zscore(list,PairwiseListMatrix([2.], false))
        @test_throws ErrorException PairwiseListMatrices.zscore(list,PairwiseListMatrix([2.,2.,2.,2.,2.,2.],true))
    end
end

@testset "Macros" begin

    PLM = PairwiseListMatrix([1,2,3], false)

    @iteratelist PLM Base.Test.@test list[k] == k

    list_values = [1,2,3,4,5,6]
    PLMtrue  = PairwiseListMatrix(list_values, true)
    PLMfalse = PairwiseListMatrix(list_values, false)
    full_t   = full(PLMtrue)
    full_f   = full(PLMfalse)

    @iteratelist PLMtrue  Base.Test.@test list[k] == :($list_values)[k]
    @iteratelist PLMfalse Base.Test.@test list[k] == :($list_values)[k]

    @iteratediag PLMtrue  Base.Test.@test false
    @iteratediag PLMfalse Base.Test.@test diag[k] == 0

    @iterateupper PLMtrue  true  list[k] = :($list_values)[k]
    @iterateupper PLMfalse false list[k] = :($list_values)[k]

    @iterateupper PLMtrue  true  list[k] = :($full_t)[i,j]
    @iterateupper PLMtrue  false list[k] = :($full_t)[i,j]

    @iterateupper PLMtrue  true  list[k] = :($full_f)[i,j]
    @iterateupper PLMfalse false list[k] = :($full_f)[i,j]
end

@testset "IO" begin

    @testset "to/from table" begin

        table = [ "A" "B" 10
                  "A" "C" 20
                  "B" "C" 30 ]

        list = setlabels(PairwiseListMatrix([10,20,30]), ["A", "B", "C"])
        list_diag = setlabels(PairwiseListMatrix([0,10,20,0,30,0], true), ["A", "B", "C"])

        @test to_table(list, diagonal=false) == table
        @test from_table(table, false, diagonalvalue=0) == list

        @test to_table(list) == [ "A" "A" 0
                                  "A" "B" 10
                                  "A" "C" 20
                                  "B" "B" 0
                                  "B" "C" 30
                                  "C" "C" 0]

        @test to_table(list) == to_table(list_diag)
        @test to_table(list, diagonal=false) == to_table(list_diag, diagonal=false)
    end

    @testset "DataFrames" begin

        values   = [1,2,3,4,5,6]
        PLMtrue  = PairwiseListMatrix(values, true)
        PLMfalse = PairwiseListMatrix(values, false)

        df = DataFrame(to_dict(PLMtrue))
        @test PLMtrue == from_table(df, true)

        df = DataFrame(to_dict(PLMtrue, diagonal=false))
        @test triu(PLMtrue,1) == triu(from_table(df, false).array,1)

        df = DataFrame(to_dict(PLMfalse))
        @test PLMfalse == from_table(df, true)

        df = DataFrame(to_dict(PLMfalse, diagonal=false))
        @test PLMfalse == from_table(df, false)
    end

    @testset "Write" begin

        values   = [1,2,3,4,5,6]
        PLMtrue  = PairwiseListMatrix(values, true)
        PLMfalse = PairwiseListMatrix(values, false)
        filename = joinpath(tempdir(), "pairwiselistmatrices_test.temp")

        try
          writecsv(filename, PLMtrue, diagonal=true)
          @test map(string,readcsv(filename,Int)) == map(string,to_table(PLMtrue,diagonal=true))
        finally
          rm(filename)
        end

        try
          writecsv(filename,PLMtrue,diagonal=false)
          @test map(string,readcsv(filename,Int)) == map(string,to_table(PLMtrue,diagonal=false))
        finally
          rm(filename)
        end

        try
          writecsv(filename,PLMfalse,diagonal=true)
          @test map(string,readcsv(filename,Int)) == map(string,to_table(PLMfalse,diagonal=true))
        finally
          rm(filename)
        end

        try
          writecsv(filename,PLMfalse,diagonal=false)
          @test map(string,readcsv(filename,Int)) == map(string,to_table(PLMfalse,diagonal=false))
        finally
          rm(filename)
        end

    end
end
