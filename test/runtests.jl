using PairwiseListMatrices
using DataFrames
using NamedArrays

using Test
using DelimitedFiles
using LinearAlgebra
using Statistics

@testset "PairwiseListMatrices" begin

    list = [1,-2,3]
    mat_diag_triu =  [ 1 -2
                        0 3 ]
    mat_diag_tril =  [ 1 0
                       -2 3 ]
    mat_triu = [ 0 1 -2
                  0 0 3
                  0 0 0 ]
    mat_tril = [ 0 0 0
                  1 0 0
                 -2 3 0 ]
    mat_diag_sym = [ 1 -2
                     -2 3 ]
    mat_sym = [ 0 1 -2
                 1 0 3
                -2 3 0 ]

    plmd_triu = UpperTriangular(PairwiseListMatrix(list, true))
    plmd_tril = LowerTriangular(PairwiseListMatrix(list, true))
    plm_triu  = UpperTriangular(PairwiseListMatrix(list))
    plm_tril  = LowerTriangular(PairwiseListMatrix(list))
    plmd_sym  = PairwiseListMatrix(list, true)
    plm_sym   = PairwiseListMatrix(list)

    @testset "Creation from list" begin

        @test plmd_triu == mat_diag_triu
        @test plmd_tril == mat_diag_tril
        @test plm_triu  == mat_triu
        @test plm_tril  == mat_tril
        @test plmd_sym  == mat_diag_sym
        @test plm_sym   == mat_sym
    end

    @testset "Matrix" begin

        @test Matrix(plmd_triu) == mat_diag_triu
        @test Matrix(plmd_tril) == mat_diag_tril
        @test Matrix(plm_triu ) == mat_triu
        @test Matrix(plm_tril ) == mat_tril
        @test Matrix(plmd_sym ) == mat_diag_sym
        @test Matrix(plm_sym  ) == mat_sym
    end

    @testset "Getters" begin

        @test getlist(plmd_sym) == plmd_sym.list
        @test getlist(plm_sym)  ==  plm_sym.list

        @test getdiag(plmd_sym) == plmd_sym.diag
        @test getdiag(plm_sym)  ==  plm_sym.diag
    end

    @testset "Length" begin

        @test lengthlist(2, true)  == lengthlist(plmd_sym)
        @test lengthlist(3, false) == lengthlist(plm_sym)

        @test length(plmd_sym) == 4
        @test length(plm_sym) == 9
    end

    @testset "Linear indexing" begin

        for i in 1:length(plm_sym)
            @test plm_sym[i] == mat_sym[i]
            plm_sym[i] = 100
            @test plm_sym[i] == 100
            plm_sym[i] = mat_sym[i]
        end
        for i in 1:length(plmd_sym)
            @test plmd_sym[i] == mat_diag_sym[i]
            plmd_sym[i] = 100
            @test plmd_sym[i] == 100
            plmd_sym[i] = mat_diag_sym[i]
        end

    end

    @testset "Convert" begin

        @test mat_diag_sym ==
                convert(PairwiseListMatrix{Int, true, Vector{Int}},  mat_diag_sym)
        @test mat_diag_sym ==
                convert(PairwiseListMatrix{Int, false, Vector{Int}}, mat_diag_sym)

        @test convert(Matrix{Float64}, mat_sym) ==
                convert(PairwiseListMatrix{Float64, true, Vector{Float64}},  mat_sym)
        @test convert(Matrix{Float64}, mat_sym) ==
                convert(PairwiseListMatrix{Float64, false, Vector{Float64}}, mat_sym)

        @test convert(Matrix{Int}, plm_sym)  == mat_sym
        @test convert(Matrix{Int}, plmd_sym) == mat_diag_sym

        @test convert(PairwiseListMatrix{Int, false, Vector{Int}}, plmd_sym) == plmd_sym
        @test convert(PairwiseListMatrix{Int, true, Vector{Int}}, plm_sym)   == plm_sym
    end

    @testset "Unary operations" begin

        for f in [ -, abs, x -> sqrt(abs(x))]
            @test broadcast(f, plmd_triu) == broadcast(f, mat_diag_triu)
            @test broadcast(f, plmd_tril) == broadcast(f, mat_diag_tril)
            @test broadcast(f, plm_triu ) == broadcast(f, mat_triu)
            @test broadcast(f, plm_tril ) == broadcast(f, mat_tril)
            @test broadcast(f, plmd_sym ) == broadcast(f, mat_diag_sym)
            @test broadcast(f, plm_sym  ) == broadcast(f, mat_sym)
        end
    end

    @testset "Binary operations" begin

        for f in [ -, +, * ]
            @test broadcast(f, plmd_triu, plmd_triu) == broadcast(f, mat_diag_triu, mat_diag_triu)
            @test broadcast(f, plmd_tril, plmd_tril) == broadcast(f, mat_diag_tril, mat_diag_tril)
            @test broadcast(f, plm_triu , plm_triu ) == broadcast(f, mat_triu,      mat_triu     )
            @test broadcast(f, plm_tril , plm_tril ) == broadcast(f, mat_tril,      mat_tril     )
            @test broadcast(f, plmd_sym , plmd_sym ) == broadcast(f, mat_diag_sym,  mat_diag_sym )
            @test broadcast(f, plm_sym  , plm_sym  ) == broadcast(f, mat_sym,       mat_sym      )
        end

        @test plmd_sym ./ plmd_sym == mat_diag_sym ./ mat_diag_sym

        # https://github.com/JuliaLang/julia/issues/19615
        for f in [ *, / ]
            @test broadcast(f, plmd_triu, 2) == broadcast(f, mat_diag_triu, 2)
            @test broadcast(f, plmd_tril, 2) == broadcast(f, mat_diag_tril, 2)
            @test broadcast(f, plm_triu , 2) == broadcast(f, mat_triu,      2)
            @test broadcast(f, plm_tril , 2) == broadcast(f, mat_tril,      2)
        end

        for f in [ -, +, *, / ]
            @test broadcast(f, plmd_sym , 2) == broadcast(f, mat_diag_sym,  2)
            @test broadcast(f, plm_sym  , 2) == broadcast(f, mat_sym,       2)
        end

        @testset "./ with Float64 & Int" begin

            plm_t = PairwiseListMatrix([.5, .4, .3], true)
            plm_f = PairwiseListMatrix([1., 1., 1.], false)
            fill!(plm_f, 1.0)
            mat_t = Matrix(plm_t)
            mat_f = Matrix(plm_f)

            for value in  (4, 4.0)
                result =  plm_t ./ value
                @test isa(result, PairwiseListMatrix{Float64,true,Vector{Float64}})
                @test result == (mat_t ./ value)

                result =  value ./ plm_t
                @test isa(result, PairwiseListMatrix{Float64,true,Vector{Float64}})
                @test result == (value ./ mat_t)

                result =  plm_f ./ value
                @test isa(result, PairwiseListMatrix{Float64,false,Vector{Float64}})
                @test result == (mat_f ./ value)

                result =  value ./ plm_f
                @test isa(result, PairwiseListMatrix{Float64,false,Vector{Float64}})
                @test result == (value ./ mat_f)
            end
        end

        # plm = PairwiseListMatrix(rand(500500), true)
        # @btime broadcast(Base.:/::typeof(/), $plm, 10.0);
        # @btime broadcast(Base.:/::typeof(/), $plm.list, 10.0);
        # @btime 0.5 .* $plm ./ 10.0;
        # @btime 0.5 .* $plm.list ./ 10.0;

    end

    @testset "Inplace Broadcast" begin

        for (p,m) in [(plmd_sym, mat_diag_sym), (plm_sym, mat_sym)]
            cp = deepcopy(p)
            cm = deepcopy(m)
            cp .+= 1.0
            cm .+= 1.0
            @test cp == cm
            cp .+= cp
            cm .+= cm
            @test cp == cm
        end
    end

    @testset "Transpose" begin

        @test transpose(plmd_triu) == plmd_tril == transpose(mat_diag_triu)
        @test transpose(plmd_tril) == plmd_triu == transpose(mat_diag_tril)
        @test transpose(plm_triu) == plm_tril == transpose(mat_triu)
        @test transpose(plm_tril) == plm_triu == transpose(mat_tril)

        @test adjoint(plmd_triu) == plmd_tril == adjoint(mat_diag_triu)
        @test adjoint(plmd_tril) == plmd_triu == adjoint(mat_diag_tril)
        @test adjoint(plm_triu) == plm_tril == adjoint(mat_triu)
        @test adjoint(plm_tril) == plm_triu == adjoint(mat_tril)

        @test transpose(plmd_sym) == plmd_sym == transpose(mat_diag_sym)
        @test transpose(plm_sym) == plm_sym == transpose(mat_sym)

        @test adjoint(plmd_sym) == plmd_sym == adjoint(mat_diag_sym)
        @test adjoint(plm_sym) == plm_sym == adjoint(mat_sym)

        @test transpose!(plm_sym) == plm_sym
        @test adjoint!(plm_sym) == plm_sym
    end

    @testset "Linear algebra" begin

        @test plmd_triu * plmd_triu == mat_diag_triu * mat_diag_triu
        @test plmd_tril * plmd_tril == mat_diag_tril * mat_diag_tril
        @test plm_triu  * plm_triu  == mat_triu      * mat_triu
        @test plm_tril  * plm_tril  == mat_tril      * mat_tril
        @test plmd_sym  * plmd_sym  == mat_diag_sym  * mat_diag_sym
        @test plm_sym   * plm_sym   == mat_sym       * mat_sym
        @test Symmetric(plm_sym) * Symmetric(plm_sym) == Symmetric(mat_sym) * Symmetric(mat_sym)

        @test plmd_triu / plmd_triu == mat_diag_triu / mat_diag_triu
        @test plmd_tril / plmd_tril == mat_diag_tril / mat_diag_tril
        @test plmd_sym  / plmd_sym  == mat_diag_sym  / mat_diag_sym
        @test plm_sym   / plm_sym   == mat_sym       / mat_sym

        for p in [:U, :S, :Vt]
            @test getproperty(svd(plmd_triu), p) ≈ getproperty(svd(mat_diag_triu), p)
            @test getproperty(svd(plmd_tril), p) ≈ getproperty(svd(mat_diag_tril), p)
            @test getproperty(svd(plm_triu ), p) ≈ getproperty(svd(mat_triu), p)
            @test getproperty(svd(plm_tril ), p) ≈ getproperty(svd(mat_tril), p)
            @test getproperty(svd(plmd_sym ), p) ≈ getproperty(svd(mat_diag_sym), p)
            @test getproperty(svd(plm_sym  ), p) ≈ getproperty(svd(mat_sym), p)
        end
    end

    @testset "Stats" begin

        for f in [ mean, std ]
            @test f(plmd_triu) == f(mat_diag_triu)
            @test f(plmd_tril) == f(mat_diag_tril)
            @test f(plm_triu ) == f(mat_triu)
            @test f(plm_tril ) == f(mat_tril)
            @test f(plmd_sym ) == f(mat_diag_sym)
            @test f(plm_sym  ) == f(mat_sym)

            @test f(plmd_sym, dims=1) == f(mat_diag_sym, dims=1)
            @test f(plm_sym , dims=1) == f(mat_sym, dims=1)
            @test f(plmd_sym, dims=2) == f(mat_diag_sym, dims=2)
            @test f(plm_sym , dims=2) == f(mat_sym, dims=2)
        end

        @test cor(plmd_sym) ≈ cor(mat_diag_sym)
        @test cor(plm_sym ) ≈ cor(mat_sym)

        # broadcast used in cor(...)
        @test plmd_sym .- mean(plmd_sym,dims=1) == mat_diag_sym .- mean(mat_diag_sym,dims=1)
        @test plmd_sym .- mean(plmd_sym,dims=2) == mat_diag_sym .- mean(mat_diag_sym,dims=2)
    end

    @testset "triu" begin

        @test triu(plmd_sym ) == triu(mat_diag_sym)
        @test triu(plm_sym  ) == triu(mat_sym)

        @test triu(plmd_sym,1) == triu(mat_diag_sym, 1)
        @test triu(plm_sym, 1) == triu(mat_sym, 1)

        @test_throws ErrorException triu!(plm_sym)
        @test_throws ErrorException triu!(plm_sym, 1)
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
            @test vec( mean_nodiag(list, dims=1) ) == [15., 20., 25.]
            @test vec( mean_nodiag(list, dims=2) ) == [15., 20., 25.]
            @test mean_nodiag(list) == 20.
        end
    end
end

@testset "Named PairwiseListMatrix" begin

    list =  [1,2,3]
    labels= ["a", "b", "c"]
    labels_diag = ["A", "B"]

    plm_diag = PairwiseListMatrix(list, true)
    plm = PairwiseListMatrix(list)
    nplm_diag = setlabels(plm_diag, labels_diag)
    nplm = setlabels(plm, labels)

    @testset "set and get labels" begin

        @test isa(nplm_diag, NamedArray)
        @test isa(nplm, NamedArray)
        @test isa(nplm_diag.array, PairwiseListMatrix)
        @test isa(nplm.array, PairwiseListMatrix)

        @test_throws AssertionError setlabels!(nplm_diag, labels)
        @test_throws AssertionError setlabels!(nplm, labels_diag)

        @test getlabels(nplm_diag) == labels_diag
        @test getlabels(nplm) == labels

        setlabels!(nplm, ["A","B","C"])
        @test getlabels(nplm) == ["A","B","C"]
        setlabels!(nplm, ["a","b","c"])
    end

    @testset "get/setindex" begin

        for i in 1:2
            for j in 2:3
                @test nplm[labels[i], labels[j]] == nplm[labels[j], labels[i]]
            end
        end

        @test nplm_diag["A", "B"] == nplm_diag["B", "A"]

        nplm_diag["A", "B"] = 20
        @test nplm_diag["B", "A"] == 20
        nplm_diag["B", "A"] = 2
        @test nplm_diag["A", "B"] == 2
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

    @testset "Named and not named arrays" begin

        @testset "Symmetric" begin

            @test issymmetric(plm)
            @test issymmetric(nplm)
            @test issymmetric(plm_diag)
            @test issymmetric(nplm_diag)
        end

        @testset "Diagonal" begin

            @test hasdiagonal(plm_diag)
            @test hasdiagonal(nplm_diag)
            @test !hasdiagonal(plm)
            @test !hasdiagonal(nplm)

            @test diagonal(plm_diag) == [1,3]
            @test diagonal(nplm_diag) == [1,3]
            @test diagonal(plm) == [0,0,0]
            @test diagonal(nplm) == [0,0,0]
        end

        @testset "eltype" begin

            @test eltype(plm) == Int
            @test eltype(nplm) == Int
            @test eltype(plm_diag) == Int
            @test eltype(nplm_diag) == Int
        end

        @testset "delegated functions" begin

            for F in (getlist, getdiag, Matrix, lengthlist, sum_nodiag, mean_nodiag, diagonal)
                    @test F(nplm) == F(plm)
                    @test F(nplm_diag) == F(plm_diag)
            end
        end

        @testset "map and broadcast" begin

            mat = Matrix(nplm)
            mat_diag = Matrix(nplm_diag)

            @test map(sqrt,plm) == map(sqrt,mat)
            @test map(sqrt,nplm) == map(sqrt,mat)
            @test map(sqrt,plm_diag) == map(sqrt,mat_diag)
            @test map(sqrt,nplm_diag) == map(sqrt,mat_diag)

            @test broadcast(sqrt,plm) == broadcast(sqrt,mat)
            @test broadcast(sqrt,nplm) == broadcast(sqrt,mat)
            @test broadcast(sqrt,plm_diag) == broadcast(sqrt,mat_diag)
            @test broadcast(sqrt,nplm_diag) == broadcast(sqrt,mat_diag)
        end
    end
end

@testset "Vector{PairwiseListMatrix}" begin

    list_samples = [ PairwiseListMatrix(rand(10), true) for i in 1:100 ]
    list_diag_samples = [ PairwiseListMatrix(rand(6), false) for i in 1:100 ]
    full_samples = Matrix{Float64}[ Matrix(mat) for mat in list_samples ]
    full_diag_samples = Matrix{Float64}[ Matrix(mat) for mat in list_diag_samples ]

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
    @test std(list_samples) ≈ sqrt.(mean(abs2.(full_samples .- mean(full_samples, dims=3)), dims=3))
    # @test std(list_diag_samples) ≈ std(full_diag_samples,3)[:,:,1] # It uses n - 1
    @test std(list_diag_samples) ≈ sqrt.(mean(abs2.(full_diag_samples .- mean(full_diag_samples, dims=3)), dims=3))

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

    @iteratelist PLM Main.@test(list[k] == k)

    list_values = [1,2,3,4,5,6]
    PLMtrue  = PairwiseListMatrix(list_values, true)
    PLMfalse = PairwiseListMatrix(list_values, false)
    full_t   = Matrix(PLMtrue)
    full_f   = Matrix(PLMfalse)

    @iteratelist PLMtrue  Main.@test(list[k] == :($list_values)[k])
    @iteratelist PLMfalse Main.@test(list[k] == :($list_values)[k])

    @iteratediag PLMtrue  Main.@test(false)
    @iteratediag PLMfalse Main.@test(diag[k] == 0)

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
          writedlm(filename, PLMtrue, diagonal=true)
          @test map(string,readdlm(filename,Int)) == map(string,to_table(PLMtrue,diagonal=true))
        finally
          rm(filename)
        end

        try
          writedlm(filename,PLMtrue,diagonal=false)
          @test map(string,readdlm(filename,Int)) == map(string,to_table(PLMtrue,diagonal=false))
        finally
          rm(filename)
        end

        try
          writedlm(filename,PLMfalse,diagonal=true)
          @test map(string,readdlm(filename,Int)) == map(string,to_table(PLMfalse,diagonal=true))
        finally
          rm(filename)
        end

        try
          writedlm(filename,PLMfalse,diagonal=false)
          @test map(string,readdlm(filename,Int)) == map(string,to_table(PLMfalse,diagonal=false))
        finally
          rm(filename)
        end

    end
end
