gap> START_TEST("Forms: matobj/quadformfields.tst");
gap> ToMatObj := m -> Matrix(IsPlistMatrixRep, DefaultFieldOfMatrix(m), m);;
gap> mat := 
> [[Z(2)^0,Z(2)^0,0*Z(2),0*Z(2)],[0*Z(2),Z(2)^0,0*Z(2),0*Z(2)], 
>  [0*Z(2),0*Z(2),0*Z(2),Z(2)^0],[0*Z(2),0*Z(2),0*Z(2),0*Z(2)]];
[ [ Z(2)^0, Z(2)^0, 0*Z(2), 0*Z(2) ], [ 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2) ], 
  [ 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ], [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ] ]
gap> form := QuadraticFormByMatrix(ToMatObj(mat));
< quadratic form >
gap> WittIndex(form);
1
gap> form := QuadraticFormByMatrix(ToMatObj(mat),GF(4));
< quadratic form >
gap> WittIndex(form);
2
gap> STOP_TEST("matobj/quadformfields.tst", 10000 );
