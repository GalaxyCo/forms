gap> START_TEST("Forms: matobj/hermitianformbymatrix.tst");
gap> ToMatObj := m -> Matrix(IsPlistMatrixRep, DefaultFieldOfMatrix(m), m);;
gap> gf := GF(3^2);
GF(3^2)
gap> mat := IdentityMat(4, gf);
[ [ Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3) ], [ 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3) ], 
  [ 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3) ], [ 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0 ] ]
gap> form := HermitianFormByMatrix( ToMatObj(mat), gf );
< hermitian form >
gap> Display(form);
Hermitian form
Gram Matrix:
 1 . . .
 . 1 . .
 . . 1 .
 . . . 1
gap> mat := [[Z(11)^0,0*Z(11),0*Z(11)],[0*Z(11),0*Z(11),Z(11)],
>     [0*Z(11),Z(11),0*Z(11)]];
[ [ Z(11)^0, 0*Z(11), 0*Z(11) ], [ 0*Z(11), 0*Z(11), Z(11) ], 
  [ 0*Z(11), Z(11), 0*Z(11) ] ]
gap> form := HermitianFormByMatrix(ToMatObj(mat),GF(121));
< hermitian form >
gap> Display(form);
Hermitian form
Gram Matrix:
  1  .  .
  .  .  2
  .  2  .
gap> STOP_TEST("matobj/hermitianformbymatrix.tst", 10000 );
