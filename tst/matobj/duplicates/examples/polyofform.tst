gap> START_TEST("Forms: matobj/polyofform.tst");
gap> ToMatObj := m -> Matrix(IsPlistMatrixRep, DefaultFieldOfMatrix(m), m);;
gap> mat := [ [ Z(8) , 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
>  [ 0*Z(2), Z(2)^0, Z(2^3)^5, 0*Z(2), 0*Z(2) ], 
>  [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ], 
>  [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), Z(2)^0 ], 
>  [ 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2) ] ];;
gap> form := QuadraticFormByMatrix(ToMatObj(mat),GF(8));
< quadratic form >
gap> PolynomialOfForm(form);
Z(2^3)*x_1^2+x_2^2+Z(2^3)^5*x_2*x_3+x_4*x_5
gap> STOP_TEST("matobj/polyofform.tst", 10000 );
