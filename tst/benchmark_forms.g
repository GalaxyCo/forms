GOMinus := function(nn, qq)
    return GeneralOrthogonalGroup(-1, nn, qq);
end;

RandomMatrixInvertible := function(n, q)
    local F, M, i, j;
    F := GF(q);
    while true do
        M := NullMat(n, n, F);
        for i in [1..n] do
            for j in [1..n] do
                M[i][j] := PseudoRandom(F);
            od;
        od;
        ConvertToMatrixRep(M, F);
        if RankMat(M) = n then
            return M;
        fi;
    od;
end;

# random invertible conjugated matrices
RandomMatsConjugated := function(n, q)
    local As, P, Pinv, Bs;
    As := List( [ 1 .. 2 ], i -> RandomMatrixInvertible(n, q) );
    P := RandomMatrixInvertible(n, q);
    Pinv := Inverse(P);
    # bs is conjugate of As (non zero solution space)
    Bs := List( As, A -> Pinv * A * P );
    return [As, Bs];
end;

#random invertible matrices
RandomMats := function(n, q)
    local As, Bs;
    As := List( [ 1 .. 2 ], i -> RandomMatrixInvertible(n, q) );
    Bs := List( [ 1 .. 2 ], i -> RandomMatrixInvertible(n, q) );
    return [As, Bs];
end;

# the generators of Sp
SpForms := function(n, q)
    local As, Bs, P;
    P := RandomMatrixInvertible(n, q);
    As := ShallowCopy(GeneratorsOfGroup(Sp(n, q)));
    As := List(As, x -> x^P);
    Bs := List(As, x -> TransposedMat(Inverse(x)));
    return [As, Bs];
end;

# random matrix and a random conjugate 
RandomSingleConjugated := function(n, q)
    local As, Bs;
    As := [RandomMatrixInvertible(n, q)];
    Bs := [As[1]^RandomMatrixInvertible(n, q)];
    return [As, Bs];
end;

# helped by claude AI
IntertwinersBenchmark := function(f)
    local sizes, n, attempts, file, results, q, F, a, As, Bs,
          t, tOwn, tMTX, resOwn, resMTX, M1, M2, r, rr;

    # settings
    n := 100;
    attempts := 5;
    file := "intertwiners_benchmark.csv";
    sizes := [ 2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 19, 23, 25, 29, 31, 37, 41, 43, 47, 53, 59, 67, 73, 81, 89, 101, 113, 125, 137, 157, 169, 191, 211, 233, 263, 289, 331, 359, 397, 439, 491, 547, 607, 673, 751, 829, 929, 1024 ];

    results := [];
    PrintTo( file, "q;own_ms;mtx_ms\n" );

    for q in sizes do
        F := GF( q );
        tOwn := 0;
        tMTX := 0;

        for a in [ 1 .. attempts ] do
            rr := f(n , q);
            As := rr[1];
            Bs := rr[2];

            t := Runtime();
            resOwn := MatrixIntertwiners( As, Bs, F );
            tOwn := tOwn + ( Runtime() - t ); 

            t := Runtime();
            M1 := GModuleByMats( As, F );
            M2 := GModuleByMats( Bs, F );
            resMTX := MTX.Homomorphisms( M1, M2 );
            tMTX := tMTX + ( Runtime() - t ); 

            if Length( resOwn ) <> Length( resMTX ) then
                Error( "dimension mismatch for q = ", q, ": ",
                       Length( resOwn ), " vs ", Length( resMTX ) );
            fi;
        od;

        r := [ q, tOwn / attempts, tMTX / attempts ];
        Add( results, r );
        AppendTo( file, r[1], ";", Int( r[2] ), ";", Int( r[3] ), "\n" );
        Print( "q = ", q, "   own: ", Int( r[2] ), " ms   MTX: ",
               Int( r[3] ), " ms\n" );
    od;

    return results;
end;


FileBenchmark := function()
    local bench_time, number_to_average, Ns, Qs, Gs, Gs_names, n, q, k, g, G, mine_total, theirs_total, i, H_conj, H, start_mine, time_after, time_start, time_after_after, average_mine, average_theirs, waste1, waste2, file_name;
    bench_time := Runtime();
    number_to_average := 1;
    Ns := [100, 300, 500];
    Qs := [5, 3^2, 2^4];
    file_name := "benchmark.csv";
    Gs := [Sp, SU, GOMinus];
    Gs_names := ["Sp", "SU", "GOMinus"];
    PrintTo(file_name, "GroupName;n;q;average_time_ours;average_time_theirs");
    AppendTo(file_name, "\n");
    for n in Ns do
        for q in Qs do
            for k in [1..Size(Gs)] do
                g := Gs[k];
                G := g(n, q);
                mine_total := 0;
                theirs_total := 0;
                for i in [1..number_to_average] do
                    H_conj := RandomMatrixInvertible(n, q);
                    H := G^H_conj;
                    start_mine := Runtime();
                    waste1 := PreservedFormsWithScalars(H);
                    time_after := Runtime() - start_mine;
                    time_start := Runtime();
                    waste2 := PreservedForms(H);
                    time_after_after := Runtime() - time_start;
                    mine_total := mine_total + time_after;
                    theirs_total := theirs_total + time_after_after;
                od;
                # '*1000' to convert to seconds
                average_mine := Float(mine_total/(number_to_average*1000));
                average_theirs := Float(theirs_total/(number_to_average*1000));
                AppendTo(file_name, Gs_names[k],";",n,";",q,";",average_mine,";",average_theirs,"\n");
            od;
        od;
    od;
    Print("Benchmark took ", Float((Runtime() - bench_time)/1000), " seconds.");
end;
