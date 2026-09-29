gap> START_TEST("Forms: matobj/istotallysingular.tst");
gap> ToMatObj := m -> Matrix(IsPlistMatrixRep, DefaultFieldOfMatrix(m), m);;
gap> mat := [[1,0,0,0,0],[0,0,0,0,1],[0,0,0,0,0],[0,0,1,0,0],[0,0,0,0,0]]*Z(8)^0;
[ [ Z(2)^0, 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
  [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ], 
  [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
  [ 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2) ], 
  [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ] ]
gap> form := QuadraticFormByMatrix(ToMatObj(mat));
< quadratic form >
gap> sub := [[Z(2)^0,0*Z(2),Z(2^3)^6,Z(2^3),Z(2^3)^3],
>        [0*Z(2),Z(2)^0,Z(2^3)^6,Z(2^3)^2,Z(2^3)]];
[ [ Z(2)^0, 0*Z(2), Z(2^3)^6, Z(2^3), Z(2^3)^3 ], 
  [ 0*Z(2), Z(2)^0, Z(2^3)^6, Z(2^3)^2, Z(2^3) ] ]
gap> IsTotallySingularSubspace(form,sub);
true
gap> STOP_TEST("matobj/istotallysingular.tst", 10000 );
