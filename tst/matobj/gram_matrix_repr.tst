#@local f, mat, form
gap> START_TEST("Forms: matobj/gram_matrix_repr.tst");

# GramMatrix must return a matrix in the same representation as the
# matrix the form was created from.

# lists in lists
gap> f := GF(5);
GF(5)
gap> mat := IdentityMat(2, f);
[ [ Z(5)^0, 0*Z(5) ], [ 0*Z(5), Z(5)^0 ] ]
gap> form := BilinearFormByMatrix(mat, f);
< bilinear form >
gap> IsMatrix(GramMatrix(form));
true
gap> form := BilinearFormByMatrix(mat);
< bilinear form >
gap> IsMatrix(GramMatrix(form));
true
gap> form := QuadraticFormByMatrix(mat, f);
< quadratic form >
gap> IsMatrix(GramMatrix(form));
true
gap> form := QuadraticFormByMatrix(mat);
< quadratic form >
gap> IsMatrix(GramMatrix(form));
true
gap> form := HermitianFormByMatrix(mat, GF(25));
< hermitian form >
gap> IsMatrix(GramMatrix(form));
true
gap> form := BilinearFormByMatrix(NullMat(2, 2, f), f);
< trivial form >
gap> IsMatrix(GramMatrix(form));
true

# matrix objects
gap> mat := Matrix(IsPlistMatrixRep, f, IdentityMat(2, f));
<2x2-matrix over GF(5)>
gap> form := BilinearFormByMatrix(mat, f);
< bilinear form >
gap> IsPlistMatrixRep(GramMatrix(form));
true
gap> form := BilinearFormByMatrix(mat);
< bilinear form >
gap> IsPlistMatrixRep(GramMatrix(form));
true
gap> form := QuadraticFormByMatrix(mat, f);
< quadratic form >
gap> IsPlistMatrixRep(GramMatrix(form));
true
gap> form := QuadraticFormByMatrix(mat);
< quadratic form >
gap> IsPlistMatrixRep(GramMatrix(form));
true
gap> form := HermitianFormByMatrix(mat, GF(25));
< hermitian form >
gap> IsPlistMatrixRep(GramMatrix(form));
true
gap> form := BilinearFormByMatrix(ZeroMatrix(2, 2, mat), f);
< trivial form >
gap> IsPlistMatrixRep(GramMatrix(form));
true
gap> STOP_TEST("matobj/gram_matrix_repr.tst", 10000);
