gap> START_TEST("Forms: matobj/discofform.tst");
gap> ToMatObj := m -> Matrix(IsPlistMatrixRep, DefaultFieldOfMatrix(m), m);;
gap> gram := InvariantQuadraticForm(GO(-1,4,5))!.matrix;;
gap> qform := QuadraticFormByMatrix(ToMatObj(gram), GF(5));
< quadratic form >
gap> DiscriminantOfForm( qform );
"nonsquare"
gap> STOP_TEST("matobj/discofform.tst", 10000 );
