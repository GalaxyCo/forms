gap> START_TEST("Forms: matobj/test_forms12.tst");
gap> ToMatObj := m -> Matrix(IsPlistMatrixRep, DefaultFieldOfMatrix(m), m);;
gap> f := GF(5);
GF(5)
gap> mat := ToMatObj(IdentityMat(3,f));
<3x3-matrix over GF(5)>
gap> form := BilinearFormByMatrix(mat,f);
< bilinear form >
gap> Display(form);
Bilinear form
Gram Matrix:
 1 . .
 . 1 .
 . . 1
gap> mat[1,2] := One(f);
Z(5)^0
gap> Display(form);
Bilinear form
Gram Matrix:
 1 . .
 . 1 .
 . . 1
gap> f := GF(5);
GF(5)
gap> mat := ToMatObj(IdentityMat(3,f));
<3x3-matrix over GF(5)>
gap> form := QuadraticFormByMatrix(mat,f);
< quadratic form >
gap> Display(form);
Quadratic form
Gram Matrix:
 1 . .
 . 1 .
 . . 1
gap> mat[1,2] := One(f);
Z(5)^0
gap> Display(form);
Quadratic form
Gram Matrix:
 1 . .
 . 1 .
 . . 1
gap> f := GF(81);
GF(3^4)
gap> mat := ToMatObj(IdentityMat(4,f));
<4x4-matrix over GF(3)>
gap> form := HermitianFormByMatrix(mat,f);
< hermitian form >
gap> Display(form);
Hermitian form
Gram Matrix:
 1 . . .
 . 1 . .
 . . 1 .
 . . . 1
gap> mat[1,2] := One(f);
Z(3)^0
gap> Display(form);
Hermitian form
Gram Matrix:
 1 . . .
 . 1 . .
 . . 1 .
 . . . 1
gap> STOP_TEST("matobj/test_forms12.tst", 10000 );
