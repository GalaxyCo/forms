gap> START_TEST("Forms: matobj/test_forms13.tst");
gap> ToMatObj := {m, F} -> Matrix(IsPlistMatrixRep, F, m);;
gap> q := 9;
9
gap> f := GF(q);
GF(3^2)
gap> dim := 3;
3
gap> v := f^dim;
( GF(3^2)^3 )
gap> mat := ToMatObj(IdentityMat(dim,f), f);
<3x3-matrix over GF(3^2)>
gap> form := QuadraticFormByMatrix(mat,f);
< quadratic form >
gap> lines := Subspaces(v,1);
Subspaces( ( GF(3^2)^3 ), 1 )
gap> matrices := List(lines,x->BasisVectors(Basis(x)));;
gap> vectors := List(matrices,x->x[1]);;
gap> results := Collected(List(vectors,x->EvaluateForm(form,Vector(x, GramMatrix(form)[1]))));; # have to convert to vector object as list times matrix obj multiplication is dubious at this point in time. See https://github.com/gap-system/gap/issues/6635 This is in a lot of the test files, namely test_forms13, test_forms14, test_forms15, test_forms16
gap> [Zero(f),(q^(dim-1)-1)/(q-1)] in results;
true
gap> results := Collected(List(matrices,x->Matrix(x, GramMatrix(form))^form));;
gap> results := List(results, x -> [Unpack(x[1]), x[2]]);; ## have to do annoying conversion, as the results are 1by1 matrices but we compare with list list 1 by 1 matrices.
gap> [[[Zero(f)]],(q^(dim-1)-1)/(q-1)] in results;
true
gap> Number(vectors,x->IsSingularVector(form,x))=(q^(dim-1)-1)/(q-1);
true
gap> Number(matrices,x->IsTotallySingularSubspace(form,Matrix(x, GramMatrix(form))))=(q^(dim-1)-1)/(q-1); # also here we have to convert to a 1by n matrix.
true
gap> STOP_TEST("matobj/test_forms13.tst", 10000 );
