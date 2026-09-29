gap> START_TEST("Forms: matobj/bg_th_ex9.tst");
gap> ToMatObj := m -> Matrix(IsPlistMatrixRep, DefaultFieldOfMatrix(m), m);;
gap> mat := [ [ Z(2^2), Z(2^2), Z(2^2), Z(2^2), Z(2^2) ], 
>    [ 0*Z(2), Z(2^2), Z(2^2)^2, 0*Z(2), Z(2)^0 ], 
>    [ 0*Z(2), 0*Z(2), Z(2)^0, Z(2)^0, Z(2)^0 ], 
>    [ 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0, Z(2)^0 ], 
>    [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ] ];;
gap> qform := QuadraticFormByMatrix(ToMatObj(mat),GF(4));
< quadratic form >
gap> IsSingularForm(qform);
false
gap> IsDegenerateForm(qform);
#I  Testing degeneracy of the *associated bilinear form*
true
gap> biform := AssociatedBilinearForm(qform);
< bilinear form >
gap> Display(biform);
Bilinear form
Gram Matrix:
z = Z(4)
   . z^1 z^1 z^1 z^1
 z^1   . z^2   .   1
 z^1 z^2   .   1   1
 z^1   .   1   .   1
 z^1   1   1   1   .
gap> IsDegenerateForm(biform);
true
gap> STOP_TEST("matobj/bg_th_ex9.tst", 10000 );
