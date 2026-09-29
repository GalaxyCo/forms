gap> START_TEST("Forms: matobj/scalarfromsim.tst");
gap> ToMatObj := m -> Matrix(IsPlistMatrixRep, DefaultFieldOfMatrix(m), m);;
gap> gram := [ [ 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ], 
>   [ Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3) ], 
>   [ 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3) ], 
>   [ 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3), 0*Z(3), 0*Z(3) ], 
>   [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0 ], 
>   [ 0*Z(3), 0*Z(3), 0*Z(3), 0*Z(3), Z(3)^0, 0*Z(3) ] ];;
gap> form := BilinearFormByMatrix( ToMatObj(gram), GF(3) );
< bilinear form >
gap> m := [ [ Z(3)^0, Z(3)^0, Z(3), 0*Z(3), Z(3)^0, Z(3) ], 
>   [ Z(3), Z(3), Z(3)^0, 0*Z(3), Z(3)^0, Z(3) ], 
>   [ 0*Z(3), Z(3), 0*Z(3), Z(3), 0*Z(3), 0*Z(3) ], 
>   [ 0*Z(3), Z(3), Z(3)^0, Z(3), Z(3), Z(3) ], 
>   [ Z(3)^0, Z(3)^0, Z(3), Z(3), Z(3)^0, Z(3)^0 ], 
>   [ Z(3)^0, 0*Z(3), Z(3), Z(3)^0, Z(3), Z(3) ] ];;
gap> ScalarOfSimilarity( m, form );
Z(3)
gap> STOP_TEST("matobj/scalarfromsim.tst", 10000 );
