#@local filt, f, dim, ToMatObj, v, w, frob, M, bil, her, quad, triv, id
gap> START_TEST("Forms: /matobj/evaluate_form_respects_matrix_obj.tst");
gap> ToMatObj := m -> Matrix(IsGenericMatrixRep, GF(9), m);;
gap> M := ToMatObj([[Z(9), Z(9)^2, Z(9)^3], [One(GF(9)), Z(9)^5, Zero(GF(9))]]);;
gap> v := RowsOfMatrix(M)[1];; w := RowsOfMatrix(M)[2];;
gap> frob := FrobeniusAutomorphism(GF(9));; id := IdentityMapping(GF(9));;
gap> bil := BilinearFormByMatrix(ToMatObj(IdentityMat(3, GF(9))), GF(9));;
gap> her := HermitianFormByMatrix(ToMatObj(IdentityMat(3, GF(9))), GF(9));;
gap> quad := QuadraticFormByMatrix(ToMatObj(IdentityMat(3, GF(9))), GF(9));;
gap> triv := BilinearFormByMatrix(ToMatObj(NullMat(3, 3, GF(9))), GF(9));;
gap> IsTrivialForm(triv);
true
gap> v^frob;
NewVector(IsPlistVectorRep,GF(9),[ Z(3^2)^3, Z(3^2)^6, Z(3^2) ])
gap> ConstructingFilter(v^frob) = ConstructingFilter(v);
true
gap> v^id;
NewVector(IsPlistVectorRep,GF(9),[ Z(3^2), Z(3^2)^2, Z(3^2)^3 ])
gap> v^id = v;
true
gap> M^frob;
NewMatrix(IsGenericMatrixRep,GF(9),3,[ [ Z(3^2)^3, Z(3^2)^6, Z(3^2) ], [ Z(3)^0, Z(3^2)^7, 0*Z(3) ] ])
gap> ConstructingFilter(M^frob) = ConstructingFilter(M);
true
gap> M^id;
NewMatrix(IsGenericMatrixRep,GF(9),3,[ [ Z(3^2), Z(3^2)^2, Z(3^2)^3 ], [ Z(3)^0, Z(3^2)^5, 0*Z(3) ] ])
gap> M^id = M;
true
gap> [v, w]^bil;
Z(3^2)^6
gap> [M, M]^bil; 
NewMatrix(IsGenericMatrixRep,GF(9),2,[ [ Z(3), Z(3^2)^6 ], [ Z(3^2)^6, Z(3^2)^7 ] ])
gap> [M, M]^id = M * TransposedMat(M);
true;
gap> [v, w]^her;
Z(3^2)^5
gap> [M, M]^her;
NewMatrix(IsGenericMatrixRep,GF(9),2,[ [ Z(3), Z(3^2)^5 ], [ Z(3^2)^7, 0*Z(3) ] ])
gap> [v, w]^triv;
0*Z(3)
gap> v^quad;
Z(3)
gap> M^quad;
NewMatrix(IsGenericMatrixRep,GF(9),2,[ [ Z(3), Z(3^2)^6 ], [ Z(3^2)^6, Z(3^2)^7 ] ])
gap> v^triv;
0*Z(3)
gap> STOP_TEST("Forms: /matobj/evaluate_form_respects_matrix_obj.tst");
