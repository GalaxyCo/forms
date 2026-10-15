## Was written with the help of claude AI
gap> START_TEST("Formspace: FORMS_IntersectVectorSpaces");

#
#  Helpers to see if the returned vector space is identical to the vector space computed. This computes 
#

# I added this because my funciton FORMS_IntersectVectorSpaces might return a matrix with one row that is just the zero vector in case of Is2BitMatrixRep or Is8BitMatrixRep, so we have to filter this out when counting linearly independent rows..
gap> FORMS_NonZeroRows := function(A)
> if Length(A) = 0 then
>   return 0;
> fi;
> if Length(A) = 1 and IsZero(A[1]) then
>   return 0;
> fi;
> return Length(A);
> end;;
gap> FORMS_TST_Rows := function(M)
>     if IsMatrixObj(M) and not IsList(M) then
>         return Unpack(M);
>     fi;
>     return List(M, ShallowCopy);
> end;;
gap> FORMS_TST_Rank := function(M)
>     if FORMS_NonZeroRows(M) = 0 then
>         return 0;
>     fi;
>     return RankMat(M);
> end;;
gap> FORMS_TST_Convert := function(rows, F, rep)
>     local M;
>     if rep = "plist" then
>         return List(rows, r -> List(r));
>     elif rep = "compressed" then
>         M := List(rows, ShallowCopy);
>         ConvertToMatrixRep(M, F);
>         return M;
>     elif rep = "matobj" then
>         return NewMatrix(IsPlistMatrixRep, F, Length(rows[1]), rows);
>     elif rep = "generic" then
>         return NewMatrix(IsGenericMatrixRep, F, Length(rows[1]), rows);
>     fi;
>     Error("unknown representation ", rep);
> end;;
gap> 1+1;;# todo also add "generic" to this, but right now this will not work because RankMat does not work for generic matrices....
gap> FORMS_TST_Reps := ["plist", "compressed", "matobj"];;
gap> FORMS_TST_CheckOne := function(A, B, F)
>     local Arows, Brows, C, rows, n, dimU, dimW, dimSum, expected;
>     Arows := FORMS_TST_Rows(A);
>     Brows := FORMS_TST_Rows(B);
>     C := FORMS_IntersectVectorSpaces(A, B);
>     rows := FORMS_TST_Rows(C);
>     n := Length(Arows[1]);
>     dimU := FORMS_TST_Rank(Arows);
>     dimW := FORMS_TST_Rank(Brows);
>     dimSum := FORMS_TST_Rank(Concatenation(Arows, Brows));
>     expected := dimU + dimW - dimSum;
>     if FORMS_NonZeroRows(rows) <> expected then
>         Print("#E wrong number of rows: ", FORMS_NonZeroRows(rows),
>               ", expected ", expected, "\n");
>         return false;
>     fi;
>     if expected = 0 then
>         return true;
>     fi;
>     if FORMS_TST_Rank(rows) <> FORMS_NonZeroRows(rows) then
>         Print("#E rows are linearly dependent\n");
>         return false;
>     fi;
>     if FORMS_TST_Rank(Concatenation(Arows, rows)) <> dimU then
>         Print("#E result not contained in span(A)\n");
>         return false;
>     fi;
>     if FORMS_TST_Rank(Concatenation(Brows, rows)) <> dimW then
>         Print("#E result not contained in span(B)\n");
>         return false;
>     fi;
>     return true;
> end;;

#  run a case with the inputs in every representation
gap> FORMS_TST_Check := function(A, B, F)
>     local rep;
>     for rep in FORMS_TST_Reps do
>         if not FORMS_TST_CheckOne(FORMS_TST_Convert(A, F, rep),
>                                   FORMS_TST_Convert(B, F, rep), F) then
>             Print("#E failed for representation \"", rep, "\"\n");
>             return false;
>         fi;
>     od;
>     return true;
> end;;

#  return the (row) dimension of the result, for each representation
gap> FORMS_TST_Dims := function(A, B, F)
>     return List(FORMS_TST_Reps, rep -> FORMS_NonZeroRows(FORMS_TST_Rows(
>         FORMS_IntersectVectorSpaces(FORMS_TST_Convert(A, F, rep),
>                                     FORMS_TST_Convert(B, F, rep)))));
> end;;

#
#  Trivial cases
#
#  identical spaces
gap> A := IdentityMat(3, GF(2));;
gap> FORMS_TST_Check(A, A, GF(2));
true

#  same space, different bases
gap> A := [[1,0,0],[0,1,0]] * Z(5)^0;;
gap> B := [[1,1,0],[1,4,0]] * Z(5)^0;;
gap> FORMS_TST_Check(A, B, GF(5));
true

#  full space with a subspace (both orders)
gap> A := IdentityMat(4, GF(7));;
gap> B := [[1,2,3,4],[0,1,6,2]] * Z(7)^0;;
gap> FORMS_TST_Check(A, B, GF(7));
true
gap> FORMS_TST_Check(B, A, GF(7));
true

#  subspace inclusion U subspace W
gap> A := [[1,1,0,0]] * Z(3)^0;;
gap> B := [[1,0,0,0],[0,1,0,0]] * Z(3)^0;;
gap> FORMS_TST_Check(A, B, GF(3));
true

#  trivial intersection: complementary subspaces, result must be empty
gap> A := [[1,0,0,0],[0,1,0,0]] * Z(2)^0;;
gap> B := [[0,0,1,0],[0,0,0,1]] * Z(2)^0;;
gap> FORMS_TST_Check(A, B, GF(2));
true
gap> FORMS_TST_Dims(A, B, GF(2));
[ 0, 0, 0 ]

#  trivial intersection, not a direct sum of the whole space
gap> A := [[1,0,0,0,0]] * Z(7)^0;;
gap> B := [[0,1,0,0,0],[0,0,1,0,0]] * Z(7)^0;;
gap> FORMS_TST_Check(A, B, GF(7));
true

#  two equal lines given by different scalar multiples
gap> A := [[1,2,3]] * Z(11)^0;;
gap> B := [[9,7,5]] * Z(11)^0;;
gap> FORMS_TST_Check(A, B, GF(11));
true
gap> FORMS_TST_Dims(A, B, GF(11));
[ 1, 1, 1 ]

#  two distinct lines
gap> A := [[1,0]] * Z(2)^0;;
gap> B := [[1,1]] * Z(2)^0;;
gap> FORMS_TST_Check(A, B, GF(2));
true

#  ambient dimension 1
gap> FORMS_TST_Check([[Z(5)]], [[Z(5)^2]], GF(5));
true

#
#  Non-trivial cases
#
#  two planes in GF(5)^3 meeting in a line
gap> A := [[1,0,0],[0,1,0]] * Z(5)^0;;
gap> B := [[1,1,1],[0,1,1]] * Z(5)^0;;
gap> FORMS_TST_Check(A, B, GF(5));
true
gap> FORMS_TST_Dims(A, B, GF(5));
[ 1, 1, 1 ]

#  two 3-spaces in GF(7)^4 meeting in a plane
gap> A := [[1,2,0,1],[0,1,1,1],[3,0,1,0]] * Z(7)^0;;
gap> B := [[1,3,1,2],[3,1,2,1],[1,1,1,1]] * Z(7)^0;;
gap> FORMS_TST_Check(A, B, GF(7));
true

#  GF(2): hyperplanes x1+x2=0 and x3+x4=0 in GF(2)^4
gap> A := [[1,1,0,0],[0,0,1,0],[0,0,0,1]] * Z(2)^0;;
gap> B := [[1,0,0,0],[0,1,0,0],[0,0,1,1]] * Z(2)^0;;
gap> FORMS_TST_Check(A, B, GF(2));
true
gap> FORMS_TST_Dims(A, B, GF(2));
[ 2, 2, 2 ]

#  non-prime field GF(9)
gap> A := [[1, Z(9), 0, 0], [0, 0, 1, Z(9)^3]] * Z(9)^0;;
gap> B := [[1, Z(9), 1, Z(9)^3], [0, 1, 0, 0]] * Z(9)^0;;
gap> FORMS_TST_Check(A, B, GF(9));
true
gap> STOP_TEST("Formspace: FORMS_IntersectVectorSpaces");
