gap> START_TEST("Forms: matobj/test_bg_th_ex2.tst");
gap> ToMatObj := m -> Matrix(IsPlistMatrixRep, DefaultFieldOfMatrix(m), m);;
gap> mat := [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1]]*Z(9)^0;
[ [ Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3) ], [ 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3) ], 
  [ 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3) ], [ 0*Z(3), 0*Z(3), 0*Z(3), Z(3) ] ]
gap> form := BilinearFormByMatrix(ToMatObj(mat),GF(9));
< bilinear form >
gap> Display(form);
Bilinear form
Gram Matrix:
 1 . . .
 . 1 . .
 . . 1 .
 . . . 2
gap> IsReflexiveForm(form);
true
gap> IsSymmetricForm(form);
true
gap> IsAlternatingForm(form);
false
gap> r := RadicalOfForm(form);;
gap> Dimension(r);
0
gap> STOP_TEST("matobj/test_bg_th_ex2.tst", 10000 );
