using PairwiseListMatrices
using Base.Test

print("""

Test PairwiseListMatrices
=======================
""")

let list = [1,-2,3],
    list_diag_triu =  [ 1 -2
                        0 3 ],
    list_diag_tril =  [ 1 0
                       -2 3 ],
    list_triu = [ 0 1 -2
                  0 0 3
                  0 0 0 ],
    list_tril = [ 0 0 0
                  1 0 0
                 -2 3 0 ],
    list_diag_sym = [ 1 -2
                     -2 3 ],
    list_sym = [ 0 1 -2
                 1 0 3
                -2 3 0 ]

  print("""

  Test creations from list
  ------------------------
  """)

  pld_triu = PairwiseListDiagonalSquareTriangular(list)
  pld_tril = PairwiseListDiagonalSquareTriangular(list, false)
  pl_triu  = PairwiseListSquareTriangular(list)
  pl_tril  = PairwiseListSquareTriangular(list, false)
  pld_sym  = PairwiseListDiagonalSymmetric(list)
  pl_sym   = PairwiseListSymmetric(list)

  @test pld_triu == list_diag_triu
  @test pld_tril == list_diag_tril
  @test pl_triu  == list_triu
  @test pl_tril  == list_tril
  @test pld_sym  == list_diag_sym
  @test pl_sym   == list_sym

  print("""

  Test unary operations
  ---------------------
  """)

  @test -pld_triu == -list_diag_triu
  @test -pld_tril == -list_diag_tril
  @test -pl_triu  == -list_triu
  @test -pl_tril  == -list_tril
  @test -pld_sym  == -list_diag_sym
  @test -pl_sym   == -list_sym

  @test abs(pld_triu) == abs(list_diag_triu)
  @test abs(pld_tril) == abs(list_diag_tril)
  @test abs(pl_triu ) == abs(list_triu)
  @test abs(pl_tril ) == abs(list_tril)
  @test abs(pld_sym ) == abs(list_diag_sym)
  @test abs(pl_sym  ) == abs(list_sym)

  print("""

  Test binary operations
  ----------------------
  """)

  @test pld_triu - pld_triu == list_diag_triu - list_diag_triu
  @test pld_tril - pld_tril == list_diag_tril - list_diag_tril
  @test pl_triu  - pl_triu  == list_triu - list_triu
  @test pl_tril  - pl_tril  == list_tril - list_tril
  @test pld_sym  - pld_sym  == list_diag_sym - list_diag_sym
  @test pl_sym   - pl_sym   == list_sym - list_sym

  @test pld_triu + pld_triu == list_diag_triu + list_diag_triu
  @test pld_tril + pld_tril == list_diag_tril + list_diag_tril
  @test pl_triu  + pl_triu  == list_triu + list_triu
  @test pl_tril  + pl_tril  == list_tril + list_tril
  @test pld_sym  + pld_sym  == list_diag_sym + list_diag_sym
  @test pl_sym   + pl_sym   == list_sym + list_sym

  @test pld_triu .* pld_triu == list_diag_triu .* list_diag_triu
  @test pld_tril .* pld_tril == list_diag_tril .* list_diag_tril
  @test pl_triu  .* pl_triu  == list_triu .* list_triu
  @test pl_tril  .* pl_tril  == list_tril .* list_tril
  @test pld_sym  .* pld_sym  == list_diag_sym .* list_diag_sym
  @test pl_sym   .* pl_sym   == list_sym .* list_sym

  print("""

  Transpose
  ---------
  """)

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

end

print("""

List with Labels
----------------
""")

let list = [1,2,3], labels=['a', 'b', 'c'], labels_diag = ['A', 'B'],
    list_diag_triu =  PairwiseListDiagonalSquareTriangular(list, true, labels_diag),
    list_diag_tril =  PairwiseListDiagonalSquareTriangular(list, false, labels_diag),
    list_triu = PairwiseListSquareTriangular(list, true, labels),
    list_tril = PairwiseListSquareTriangular(list, false, labels),
    list_diag_sym = PairwiseListDiagonalSymmetric(list, labels_diag),
    list_sym = PairwiseListSymmetric(list, labels)

  print("""
  getlabel & getindex
  """)

  @test getlabel(list_diag_triu, 'A', 'B') == list_diag_triu[1,2] == 2
  @test getlabel(list_diag_triu, 'B', 'A') == list_diag_triu[2,1] == 0

  @test getlabel(list_diag_tril, 'A', 'B') == list_diag_tril[1,2] == 0
  @test getlabel(list_diag_tril, 'B', 'A') == list_diag_tril[2,1] == 2

  @test getlabel(list_triu, 'a', 'b') == list_triu[1,2] == 1
  @test getlabel(list_triu, 'b', 'a') == list_triu[2,1] == 0

  @test getlabel(list_tril, 'a', 'b') == list_tril[1,2] == 0
  @test getlabel(list_tril, 'b', 'a') == list_tril[2,1] == 1

  @test getlabel(list_diag_sym, 'A', 'B') == list_diag_sym[1,2] == 2
  @test getlabel(list_diag_sym, 'B', 'A') == list_diag_sym[2,1] == 2

  @test getlabel(list_sym, 'a', 'b') == list_sym[1,2] == 1
  @test getlabel(list_sym, 'b', 'a') == list_sym[2,1] == 1

  print("""
  setlabel! & setindex!
  """)

  @test_throws BoundsError setlabel!(list_diag_triu, 10, 'B', 'A')
  @test_throws BoundsError setlabel!(list_diag_tril, 10, 'A', 'B')
  @test_throws BoundsError setlabel!(list_triu, 10, 'b', 'a')
  @test_throws BoundsError setlabel!(list_tril, 10, 'a', 'b')

  setlabel!(list_diag_triu, 10, 'A', 'B')
  setlabel!(list_diag_tril, 10, 'B', 'A')
  setlabel!(list_triu, 10, 'a', 'b')
  setlabel!(list_tril, 10, 'b', 'a')

  @test list_diag_triu[1,2] == 10
  @test list_diag_tril[2,1] == 10
  @test list_triu[1,2] == 10
  @test list_tril[2,1] == 10

  list_diag_sym[1,2] = 10
  @test list_diag_sym[2,1] == 10

  list_sym[1,2] = 10
  @test list_sym[2,1] == 10

end
