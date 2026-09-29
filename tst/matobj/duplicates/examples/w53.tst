gap> START_TEST("Forms: matobj/w53.tst");
gap> ToMatObj := m -> Matrix(IsPlistMatrixRep, DefaultFieldOfMatrix(m), m);;
gap> f := GF(3);
GF(3)
gap> gram := [
> [0,0,0,1,0,0], 
> [0,0,0,0,1,0],
> [0,0,0,0,0,1],
> [-1,0,0,0,0,0],
> [0,-1,0,0,0,0],
> [0,0,-1,0,0,0]] * One(f);;
gap> form := BilinearFormByMatrix( ToMatObj(gram), f );
< bilinear form >
gap> IsSymplecticForm( form );
true
gap> Display( form );
Symplectic form
Gram Matrix:
 . . . 1 . .
 . . . . 1 .
 . . . . . 1
 2 . . . . .
 . 2 . . . .
 . . 2 . . .
gap> b := BaseChangeToCanonical( form );
< immutable compressed matrix 6x6 over GF(3) >
gap> Display( b );
 1 . . . . .
 . . . 1 . .
 . . 1 . . .
 . . . . . 1
 . . . . 1 .
 . 2 . . . .
gap> Display( b * gram * TransposedMat(b) );
 . 1 . . . .
 2 . . . . .
 . . . 1 . .
 . . 2 . . .
 . . . . . 1
 . . . . 2 .
gap> STOP_TEST("matobj/w53.tst", 10000 );
