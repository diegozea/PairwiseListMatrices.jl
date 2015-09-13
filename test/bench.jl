using PairedListMatrices
using Benchmarks
using Base.Tests

# Similar to pairwise of Distances.jl
function distances_pairwise(vecs)
  n = length(vecs)
  mat = Array(Float64, n, n)
  @inbounds for i in 1:n
    vec_i = vecs[i]
    for j in i:n
      mat[i,j] = cor(vec_i, vecs[j])
    end
    for j = 1 : (i-1)
      @inbounds mat[i,j] = mat[j,i]  # leveraging the symmetry
    end
  end
  mat
end

function using_full_symmetric(vecs)
  n = length(vecs)
  mat = Array(Float64, n, n)
  @inbounds for i in 1:n
    vec_i = vecs[i]
    for j in i:n
      mat[i,j] = cor(vec_i, vecs[j])
    end
  end
  Symmetric(mat)
end

function using_full_symmetric(vecs)
  n = length(vecs)
  mat = Array(Float64, n, n)
  @inbounds for i in 1:n
    vec_i = vecs[i]
    for j in i:n
      mat[i,j] = cor(vec_i, vecs[j])
    end
  end
  Symmetric(mat)
end

function using_sparse_symmetric(vecs)
  n = length(vecs)
  mat = spzeros(Float64, n, n)
  @inbounds for i in 1:n
    vec_i = vecs[i]
    for j in i:n
      mat[i,j] = cor(vec_i, vecs[j])
    end
  end
  Symmetric(mat)
end

function pairedlist_matindex(vecs)
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

function pairedlist_listindex(vecs)
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

const testset_small=[ rand(3) for n in 1:100 ];

@test all(using_sparse_symmetric(testset_small) .== using_full_symmetric(testset_small))
@test all(using_full_symmetric(testset_small) .== pairedlist_matindex(testset_small))
@test all(pairedlist_matindex(testset_small) .== pairedlist_listindex(testset_small))
@test all(pairedlist_listindex(testset_small) .== distances_pairwise(testset_small))

@benchmark distances_pairwise(testset_small)
@benchmark using_full_symmetric(testset_small)
@benchmark using_sparse_symmetric(testset_small)
@benchmark pairedlist_matindex(testset_small)
@benchmark pairedlist_listindex(testset_small)

const testset_big=[ rand(3) for n in  1:1000 ];

@test all(using_sparse_symmetric(testset_big) .== using_full_symmetric(testset_big))
@test all(using_full_symmetric(testset_big) .== pairedlist_matindex(testset_big))
@test all(pairedlist_matindex(testset_big) .== pairedlist_listindex(testset_big))
@test all(pairedlist_listindex(testset_big) .== distances_pairwise(testset_big))

@benchmark distances_pairwise(testset_big)
@benchmark using_full_symmetric(testset_big)
@benchmark using_sparse(testset_big)
@benchmark pairedlist_matindex(testset_big)
@benchmark pairedlist_listindex(testset_big)
