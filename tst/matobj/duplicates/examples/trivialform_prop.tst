gap> START_TEST("Forms: matobj/trivialform_prop.tst");
gap> ToMatObj := m -> Matrix(IsPlistMatrixRep, DefaultFieldOfMatrix(m), m);;
gap> mat := [[0,0,0],[0,0,0],[0,0,0]]*Z(11)^0;
[ [ 0*Z(11), 0*Z(11), 0*Z(11) ], [ 0*Z(11), 0*Z(11), 0*Z(11) ], 
  [ 0*Z(11), 0*Z(11), 0*Z(11) ] ]
gap> form := QuadraticFormByMatrix(ToMatObj(mat),GF(121));
< trivial form >
gap> IsReflexiveForm(form);
true
gap> IsAlternatingForm(form);
true
gap> IsSymmetricForm(form);
true
gap> IsOrthogonalForm(form);
false
gap> IsPseudoForm(form);
false
gap> IsSymplecticForm(form);
true
gap> IsDegenerateForm(form);
true
gap> IsSingularForm(form);
true
gap> BaseField(form);
GF(11^2)
gap> GramMatrix(form);
[ [ 0*Z(11), 0*Z(11), 0*Z(11) ], [ 0*Z(11), 0*Z(11), 0*Z(11) ], 
  [ 0*Z(11), 0*Z(11), 0*Z(11) ] ]
gap> RadicalOfForm(form);
<vector space of dimension 3 over GF(11^2)>
gap> STOP_TEST("matobj/trivialform_prop.tst", 10000 );
