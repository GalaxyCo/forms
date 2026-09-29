gap> START_TEST("Forms: matobj/bg_th_ex4.tst");
gap> ToMatObj := m -> Matrix(IsPlistMatrixRep, DefaultFieldOfMatrix(m), m);;
gap> mat := [[0,1,0,0,0,0],[1,0,0,0,0,0],[0,0,0,0,0,1],
>         [0,0,0,0,1,0],[0,0,0,1,0,0],[0,0,1,0,0,0]]*Z(16)^0;
[ [ 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
  [ Z(2)^0, 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
  [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ], 
  [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2) ], 
  [ 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2) ], 
  [ 0*Z(2), 0*Z(2), Z(2)^0, 0*Z(2), 0*Z(2), 0*Z(2) ] ]
gap> form := BilinearFormByMatrix(ToMatObj(mat),GF(16));
< bilinear form >
gap> Display(form);
Bilinear form
Gram Matrix:
 . 1 . . . .
 1 . . . . .
 . . . . . 1
 . . . . 1 .
 . . . 1 . .
 . . 1 . . .
gap> IsSymmetricForm(form);
true
gap> IsAlternatingForm(form);
true
gap> IsDegenerateForm(form);
false
gap> WittIndex(form);
3
gap> STOP_TEST("matobj/bg_th_ex4.tst", 10000 );
