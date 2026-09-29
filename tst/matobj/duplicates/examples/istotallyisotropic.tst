gap> START_TEST("Forms: matobj/istotallyisotropic.tst");
gap> ToMatObj := m -> Matrix(IsPlistMatrixRep, DefaultFieldOfMatrix(m), m);;
gap> mat := [[1,0,0,0],[0,-1,0,0],[0,0,0,1],[0,0,1,0]]*Z(7)^0;
[ [ Z(7)^0, 0*Z(7), 0*Z(7), 0*Z(7) ], [ 0*Z(7), Z(7)^3, 0*Z(7), 0*Z(7) ], 
  [ 0*Z(7), 0*Z(7), 0*Z(7), Z(7)^0 ], [ 0*Z(7), 0*Z(7), Z(7)^0, 0*Z(7) ] ]
gap> form := BilinearFormByMatrix(ToMatObj(mat));
< bilinear form >
gap> sub:= [[Z(7)^0,0*Z(7),Z(7)^0,Z(7)],[0*Z(7),Z(7)^0,Z(7)^0,Z(7)^4]];
[ [ Z(7)^0, 0*Z(7), Z(7)^0, Z(7) ], [ 0*Z(7), Z(7)^0, Z(7)^0, Z(7)^4 ] ]
gap> IsTotallyIsotropicSubspace(form,sub);
true
gap> mat := IdentityMat(6,GF(2));
[ <a GF2 vector of length 6>, <a GF2 vector of length 6>, 
  <a GF2 vector of length 6>, <a GF2 vector of length 6>, 
  <a GF2 vector of length 6>, <a GF2 vector of length 6> ]
gap> form := HermitianFormByMatrix(ToMatObj(mat),GF(4));
< hermitian form >
gap> sub := [[Z(2)^0,0*Z(2),0*Z(2),Z(2)^0,Z(2)^0,Z(2)^0], 
>   [0*Z(2),Z(2)^0,0*Z(2),Z(2^2)^2,Z(2^2),Z(2)^0], 
>   [0*Z(2),0*Z(2),Z(2)^0,Z(2)^0,Z(2^2),Z(2^2)^2]];
[ [ Z(2)^0, 0*Z(2), 0*Z(2), Z(2)^0, Z(2)^0, Z(2)^0 ], 
  [ 0*Z(2), Z(2)^0, 0*Z(2), Z(2^2)^2, Z(2^2), Z(2)^0 ], 
  [ 0*Z(2), 0*Z(2), Z(2)^0, Z(2)^0, Z(2^2), Z(2^2)^2 ] ]
gap> IsTotallyIsotropicSubspace(form,sub);
true
gap> STOP_TEST("matobj/istotallyisotropic.tst", 10000 );
