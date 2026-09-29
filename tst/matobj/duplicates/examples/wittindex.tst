gap> START_TEST("Forms: matobj/wittindex.tst");
gap> ToMatObj := m -> Matrix(IsPlistMatrixRep, DefaultFieldOfMatrix(m), m);;
gap> mat := [[0,0,1,0,0],[0,0,0,0,0],[-1,0,0,0,0],[0,0,0,0,0],[0,0,0,0,0]]*Z(7)^0;
[ [ 0*Z(7), 0*Z(7), Z(7)^0, 0*Z(7), 0*Z(7) ], 
  [ 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7) ], 
  [ Z(7)^3, 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7) ], 
  [ 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7) ], 
  [ 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7) ] ]
gap> form := BilinearFormByMatrix(ToMatObj(mat),GF(7));
< bilinear form >
gap> WittIndex(form);
1
gap> RadicalOfForm(form);
<vector space of dimension 3 over GF(7)>
gap> Dimension(last);
3
gap> mat := IdentityMat(6,GF(5));
[ [ Z(5)^0, 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5) ], 
  [ 0*Z(5), Z(5)^0, 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5) ], 
  [ 0*Z(5), 0*Z(5), Z(5)^0, 0*Z(5), 0*Z(5), 0*Z(5) ], 
  [ 0*Z(5), 0*Z(5), 0*Z(5), Z(5)^0, 0*Z(5), 0*Z(5) ], 
  [ 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), Z(5)^0, 0*Z(5) ], 
  [ 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), Z(5)^0 ] ]
gap> form := QuadraticFormByMatrix(ToMatObj(mat),GF(5));
< quadratic form >
gap> WittIndex(form);
3
gap> mat := IdentityMat(6,GF(7));
[ [ Z(7)^0, 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7) ], 
  [ 0*Z(7), Z(7)^0, 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7) ], 
  [ 0*Z(7), 0*Z(7), Z(7)^0, 0*Z(7), 0*Z(7), 0*Z(7) ], 
  [ 0*Z(7), 0*Z(7), 0*Z(7), Z(7)^0, 0*Z(7), 0*Z(7) ], 
  [ 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7), Z(7)^0, 0*Z(7) ], 
  [ 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7), 0*Z(7), Z(7)^0 ] ]
gap> form := QuadraticFormByMatrix(ToMatObj(mat),GF(7));
< quadratic form >
gap> WittIndex(form);
2
gap> STOP_TEST("matobj/wittindex.tst", 10000 );
