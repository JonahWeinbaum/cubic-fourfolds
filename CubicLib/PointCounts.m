//////////////////////////////////////////////
// PointCounts.m
//
// Contains the functions necessary to run the specialized pointcounting
// code of CubicLib.

//////////////////////////
// CONSTANTS

PATH_TO_LIB := PathToLib();

k := FiniteField(2);
F4 := FiniteField(4);
basF4 := Basis(F4);
k6 := RSpace(k, 6);
k4 := RSpace(k, 4);

G_4 := GL(6, F4);
P5_4<[x]> := ProjectiveSpace(F4, 5);
R_4 := CoordinateRing(P5_4);
V_4, Bit_4 := GModule(G_4, R_4, 3);


lines := ReadLinesIndex();    
echforms := lines;

ptsonlines := AssociativeArray();
monoevaluated := AssociativeArray();

//Loop over all lines
for form in echforms do 
    nullsp := NullspaceOfTranspose(form);
    bas := Basis(nullsp);
    nullsp := ExtendField(nullsp, F4);
    pt1 := bas[1];
    pt2 := bas[2];
    pt3 := (nullsp!pt1)*1 + (nullsp!pt2)*F4.1;

    ptsonlines[form] := <pt1, pt2, pt3>;
    
    monoevaluated[form] := <0,0,0,0>;

    for mono in Basis(V_4) do 

        eval1 := Evaluate(mono @@ Bit_4, Eltseq(pt1));
        eval2 := Evaluate(mono @@ Bit_4, Eltseq(pt2));
        eval3 := Evaluate(mono @@ Bit_4, Eltseq(pt3));

        eval3bas1 := Eltseq(eval3)[1];
        eval3bas2 := Eltseq(eval3)[2];

        evals := [Integers() ! x : x in [eval1, eval2, eval3bas1, eval3bas2]];
        for i in [1..4] do
            monshiftLeft := ShiftLeft(monoevaluated[form][i], 1);
            monoevaluated[form][i] := BitwiseOr(monshiftLeft, evals[i]);
        end for;
    end for;
end for;

function paritycalc(bitstr)
    parity := false;
    while bitstr ne 0 do
        parity := not parity;
        bitstr := BitwiseAnd(bitstr, (bitstr - 1));
    end while;  
    return parity;
end function;


intrinsic LinesThrough(f :: RngMPolElt) -> SeqEnum
{Given a cubic polynomial f, return the set of lines through f.}

    f := CubicToInt(f);    
    
    lines := [];
    for form in Keys(monoevaluated) do 
        evals := monoevaluated[form];
	cont_flag := false; 

	for j in [1..4] do
	    evalj := BitwiseAnd(monoevaluated[form][j], f);
	    if paritycalc(evalj) then
		cont_flag := true;
		break;
	    end if;
	end for;
	if cont_flag then continue; end if;

        lines cat:= [form];
    end for;

    return lines;
end intrinsic;

intrinsic ConicFibrationCoefficients(cubic) -> SeqEnum
{Given a cubic containing the line x1 = x2 = x3 = x4 = 0, return the coefficients
of the conic fibration as polynomials in a new set of variables.}

    R<y0, y1, y2, y3> := PolynomialRing(k, 4);
    RR<y4, y5> := PolynomialRing(R, 4);
 
    f := Evaluate(cubic, [y0, y1, y2, y3, y4, y5]);
    A := MonomialCoefficient(f, y4^2);
    B := MonomialCoefficient(f, y4*y5);
    C := MonomialCoefficient(f, y5^2);
    D := MonomialCoefficient(f, y4);
    E := MonomialCoefficient(f, y5);
    F := MonomialCoefficient(f, 1);
    
    return [A, B, C, D, E, F];
end intrinsic;

intrinsic CubicFromConicCoefficients(conicCoeffs, R) -> RngMPol
{Produce a cubic with parent R from the conic fibration coefficients.}

    vvec := [R.i : i in [1..4]];
    x5 := R.5;
    x6 := R.6;
    mons := [x5^2, x5*x6, x6^2, x5, x6, 1];
    
    return &+[mons[i] * Evaluate(conicCoeffs[i], vvec) : i in [1..6]];    
end intrinsic;

intrinsic CorankTwoLocus(conicCoeffs) -> Sch
{Return the corank 2 locus of the conic fibration.}

    B := conicCoeffs[2];
    D := conicCoeffs[4];
    E := conicCoeffs[5];

    R := Parent(B);
    k := BaseRing(R);
    
    return Scheme(Proj(R), [B, D, E]);    
end intrinsic;

intrinsic StandardizeConicCoefficients(conicCoeffs) -> SeqEnum
{Given a list of conic fibration coefficients, try to move a singular
point of the discriminant to (0:0:0:1)}

    B := conicCoeffs[2];
    D := conicCoeffs[4];
    E := conicCoeffs[5];

    R := Parent(B);
    k := BaseRing(R);
    
    somesing := Scheme(Proj(R), [B, D, E]);
    somesingpts := Points(somesing);
        
    error if #somesingpts eq 0, "No rational points in rank 2 locus.";

    // Check that the discriminantal quintic in P3 has singular point where by B=D=E=0, 
    // and change coordinates so one of these singular points is in position (0:0:0:1).

    V4 := VectorSpace(k, 4);
    v  := V4 ! Eltseq(somesingpts[1]);
    extBase := ExtendBasis([v], V4);
    revextBase := Reverse(extBase); // Puts v at the bottom.

    M4 := Matrix(revextBase)^(-1);
    invM4 := M4^(-1);
    
    // Change coordinates to move the singular point.
    yvec := Vector(R, [R.i : i in [1..4]]);
    newcoeffs := [Evaluate(co, Eltseq(yvec * ChangeRing(invM4, R))) : co in conicCoeffs];

    return newcoeffs;
end intrinsic;


intrinsic ConicFibrationForm(cubic, line::ModMatFldElt) -> RngMPol, SeqEnum
{Given a cubic f and a line L contained in that cubic, change coordinates so that
L is the standard line y0 = y1 = y2 = y3 = 0. Also return the coefficients of the
conic fibration.}

    k := GF(2);
    V6 := VectorSpace(k, 6);
    R<y0, y1, y2, y3> := PolynomialRing(k, 4);
    RR<y4, y5> := PolynomialRing(R, 4);
    
    // Compute a transformation M such that cubic^M has a line given by y0=...=y3=0.
    vectorizedline := [V6 ! line[i] : i in [1..4]];
    M := Matrix(ExtendBasis(vectorizedline, V6))^(-1);

    // Update the cubic and extract coefficients.
    g := cubic^M;
    f := Evaluate(g, [y0, y1, y2, y3, y4, y5]);
    A := MonomialCoefficient(f, y4^2);
    B := MonomialCoefficient(f, y4*y5);
    C := MonomialCoefficient(f, y5^2);
    D := MonomialCoefficient(f, y4);
    E := MonomialCoefficient(f, y5);
    F := MonomialCoefficient(f, 1);

    return g, [A,B,C,D,E,F];
end intrinsic;

intrinsic MatrixToLine(P::Prj, mat::ModMatFldElt) -> Sch
{Converst a 4x6 matrix into a line.}
    eqns := [&+[P.i * row[i] : i in [1..6]] : row in Rows(mat)];
    L := Scheme(P, eqns);
    return L;
end intrinsic;

    
intrinsic LineToMatrix(line::Sch) -> ModFldMatElt
{}
    // First convert the line into Echelon form. We do so by extracting the
    // coefficients of the defining equations.
    assert Dimension(line) eq 1 and Degree(line) eq 1;
    eqns := DefiningEquations(line);

    R := Parent(eqns[1]);
    mat := Matrix(BaseRing(R), [[MonomialCoefficient(ll, R.i) : i in [1..6]] : ll in eqns]);
    return mat;
end intrinsic;

intrinsic ConicFibrationForm(cubic, line::Sch) -> RngMPol, SeqEnum
{}
    return ConicFibrationForm(cubic, LineToMatrix(line));
end intrinsic;


intrinsic AddingtonAuelStandardForm(cubic : Nonflat:=false) -> RngMPol, SeqEnum
{Transforms a cubic into a form as outlined in Addington-Auel Section 3.}
    
    k := GF(2);
    V6 := VectorSpace(k, 6);
    R<y0, y1, y2, y3> := PolynomialRing(k, 4);
    RR<y4, y5> := PolynomialRing(R, 4);
    
    flines := LinesThrough(cubic);

    for line in flines do
        // Compute a transformation M such that cubic^M has a line given by y0=...=y3=0.
        vectorizedline := [V6 ! line[i] : i in [1..4]];
        M := Matrix(ExtendBasis(vectorizedline, V6))^(-1);

        // Update the cubic and extract coefficients.
	g := cubic^M;
	f := Evaluate(g, [y0, y1, y2, y3, y4, y5]);
	A := MonomialCoefficient(f, y4^2);
	B := MonomialCoefficient(f, y4*y5);
	C := MonomialCoefficient(f, y5^2);
	D := MonomialCoefficient(f, y4);
	E := MonomialCoefficient(f, y5);
	F := MonomialCoefficient(f, 1);

        // Ensure that the choice is good.
	I := Saturation(ideal<R | [A, B, C, D, E, F]>) ;
	somesing := Scheme(ProjectiveSpace(R), [B, D, E]);
	somesingpts := Points(somesing);
        
        if Nonflat and #somesingpts ne 0 then
            break;
        elif 1 in I and #somesingpts ne 0 then
            break;
	end if;
    end for;

    if (1 notin I and not Nonflat) or (#somesingpts eq 0) then
	error "Cannot convert cubic into Addington-Auel standard form.";
    end if;

    //Check that the discriminantal quintic in P3 has singular point where by B=D=E=0, 
    //and change coordinates so one of these singular points is in position (0:0:0:1).

    V4 := VectorSpace(k, 4);
    v := V4 ! Eltseq(somesingpts[1]);
    extBase := ExtendBasis([v], V4);
    revextBase := Reverse(extBase); // Puts v at the bottom.

    M4 := Matrix(revextBase)^(-1);
    invM4 := M4^(-1);
    
    // Change coordinates to move the singular point.
    coeffs := [A, B, C, D, E, F];
    yvec := Vector(R, [y0, y1, y2, y3]);
    newcoeffs := [Evaluate(co, Eltseq(yvec * ChangeRing(invM4, R))) : co in coeffs];

    // Also apply the change of coordinates to g. For consistency.
    gPar  := Parent(g);
    xvec4 := Vector(gPar, [gPar.i : i in [1..4]]);

    newg := Evaluate(g, Eltseq(xvec4 * ChangeRing(invM4, gPar)) cat [gPar.5, gPar.6]);
    
    return newg, newcoeffs;
end intrinsic;

intrinsic ConicFibrationForCpp(cubic) -> SeqEnum
{Determine a conic fibration of the cubic, obtained via projecting from a line,
which is particularly amenable to input into C++. We require a the conic fibration 
of the form
    
    Ax^2 + Bxy + Cy^2 + Dx + Ey + F so that

such that the projection map (A : B : C) is the projection away from the distinguished 
singular point on the discriminant.}

    R<[x]> := Parent(cubic);
    k := BaseRing(R);

    g, conicCoeffs := AddingtonAuelStandardForm(cubic : Nonflat);

    // Find a point where A(p) = B(p) = C(p) = 0.
    P3 := Proj(Parent(conicCoeffs[1]));
    pt := Scheme(P3, conicCoeffs[1..3]);

    // assert Dimension(pt) eq 0; // TODO: Figure out code for the rank degenerate case.

    // If we are already in good coordinates, there is nothing to do.
    pt := Random(RationalPoints(pt));
    //if Eltseq(pt) eq [0,0,0,1] then
    //    return g, conicCoeffs;
    //end if;
    
    // Change coordinates so that pt = [0,0,0,1].
    trafo := Translation(P3, P3 ! [0,0,0,1], P3 ! pt);
    ABC := [Evaluate(Q, DefiningEquations(trafo)) : Q in conicCoeffs[1..3]];
    
    // Change coordinates so that A = y0, B = y1, C = y2.
    mons1 := [P3.i : i in [1..4]];
    mat := Matrix([[MonomialCoefficient(ll, m) : m in mons1] : ll in ABC] cat [[0,0,0,1]]);

    E, U  := EchelonForm(Transpose(mat));
    trafo := Automorphism(P3, U) * trafo;
    
    // Lift the transform to the cubic.
    trafoLift := [Evaluate(eqn, x[1..4]) : eqn in DefiningEquations(trafo)] cat x[5..6];

    // Compute new data.
    newCoeffs := [Evaluate(Q, DefiningEquations(trafo)) : Q in conicCoeffs];    
    newg := Evaluate(g, trafoLift);

    // Check the transform.
    // print conicCoeffs[1..3], newCoeffs[1..3];
    assert Seqset(newCoeffs[1..3]) subset {0, P3.1, P3.2, P3.3};

    // TODO: Need to force B=y1 in rank degenerate cases (if this is possible.)
    assert newCoeffs[2] eq P3.2;
    
    // Return.
    return newg, newCoeffs;
end intrinsic;


intrinsic ConicFibrationDiscriminant(conicCoeffs) -> RngMPol
{}
    A, B, C, D, E, F := Explode(conicCoeffs);
    disc := C*D^2 + A*E^2 + F*B^2 + B*D*E;
    return disc;
end intrinsic;

intrinsic DiscriminantProjectionEquations(conicCoeffs) -> SeqEnum
{}
    A, B, C, D, E, F := Explode(conicCoeffs);

    rrr<y0, y1, y2> := PolynomialRing(k, 3);
    RRR<y3> := PolynomialRing(rrr, 1);
    disc := Evaluate(C*D^2 + A*E^2 + F*B^2 + B*D*E, [RRR!y0, RRR!y1, RRR!y2, RRR!y3]);

    a := MonomialCoefficient(disc, y3^4);
    b := MonomialCoefficient(disc, y3^3);
    c := MonomialCoefficient(disc, y3^2);
    d := MonomialCoefficient(disc, y3);    
    e := MonomialCoefficient(disc, 1);
    
    return [a, b, c, d, e];
end intrinsic;


function WriteHeaderFile(fname, conicCoeffs, discCoeffs)

    A, B, C, D, E, F := Explode(conicCoeffs);
    a, b, c, d, e := Explode(discCoeffs);
    
    string :=  "#define abcde \\" cat "\n" cat 
               "    a = " cat CppInputString(a) cat ", \\" cat "\n" cat
	       "    b = " cat CppInputString(b) cat ", \\" cat "\n" cat
	       "    c = " cat CppInputString(c) cat ", \\" cat "\n" cat
	       "    d = " cat CppInputString(d) cat ", \\" cat "\n" cat
               "    e = " cat CppInputString(e) cat "\n" cat "\n" cat
	       "#define ABCDEF \\" cat "\n" cat
	       "    A = " cat CppInputString(A) cat ", \\" cat "\n" cat
	       "    B = " cat CppInputString(B) cat ", \\" cat "\n" cat
	       "    C = " cat CppInputString(C) cat ", \\" cat "\n" cat
	       "    D = " cat CppInputString(D) cat ", \\" cat "\n" cat
	       "    E = " cat CppInputString(E) cat ", \\" cat "\n" cat
	       "    F = " cat CppInputString(F) cat "\n";

    SetColumns(0);
    Write(fname, string : Overwrite);
    return 0;
end function;

intrinsic FibrationBaseSchemeOnLine(conicCoeffs) -> Sch
{}
    R := Parent(conicCoeffs[1]);
    P3 := Proj(R);
    return Scheme(P3, conicCoeffs[1..3]);
end intrinsic;

intrinsic CppCountCorrection(cppCounts, conicCoeffs) -> SeqEnum
{Add correction terms, amounting to changing P1 x P2 to the actual
exceptional fibre for the conic fibration.}

    // Determine the base scheme of the 3-dimensional conic system. This will
    // determine the correct exceptional fibres.
    Base := FibrationBaseSchemeOnLine(conicCoeffs);
    base_dim := Dimension(Base);
    
    Prng<Q> := PolynomialRing(Integers());
    skipOdd := false;

    A, B, C, D, E, F := Explode(conicCoeffs);
    
    if base_dim eq 0 then
        correctionTerm := Zero(Prng);
        
    elif base_dim eq 1 and (A eq 0 or C eq 0) then
        correctionTerm := - Q^3;

    elif base_dim eq 1 and (B eq 0) then
        correctionTerm := Zero(Prng);

    elif base_dim eq 1 then
        // Not sure if this case is correct.
        correctionTerm := Zero(Prng);
        
    elif base_dim eq 2 and ((A eq 0 and B eq 0) or (B eq 0 and C eq 0)) then
        correctionTerm := - Q^3;
        
    elif base_dim eq 2 and (A eq 0 or C eq 0) then
        correctionTerm := -2 * Q^3;

    elif base_dim eq 2 and (B eq 0) then
        correctionTerm := - Q^3;
        
    elif base_dim eq 2 and (A eq B and B eq C) then
        // In this case, the base points are quadratic conjugates.
        correctionTerm := -2 * Q^3;
        skipOdd := true;
        
    elif base_dim eq 3 then
        correctionTerm := - (Q+1) * Q;
        
    else
        error "Case not accounted for.", conicCoeffs;
    end if;

    // Odd
    if not skipOdd then
        for i in [1..#cppCounts by 2] do
            q := 2^i;
            cppCounts[i] +:= Evaluate(correctionTerm, q);        
        end for;
    end if;

    // Even
    for i in [2..#cppCounts by 2] do
        q := 2^i;
        cppCounts[i] +:= Evaluate(correctionTerm, q);        
    end for;

    return cppCounts;
end intrinsic;


intrinsic PointCounts(cubic : ExecNum:=false, Minq:=1, Maxq:=11,
                              CompilerOptimization:=false,
                              CompileEachQ:=false,                              
                              UseCache:=true,
                              PrintTimes:=false) -> SeqEnum
{Compute the pointcounts on the given cubic over Fq, where q = 2,4,...,2048.}

    // This function is mostly based off of 
    // Section 3 of Addington-Auel. See loc. cit. for details.
    //
    // The actual driver of this function is the C++ code in point_counting_cpp.

    // Now compute coefficients for the conic bundle fibration
    // that will be fed into C++ point counting code.
    if not CompileEachQ then
        g, conicCoeffs := AddingtonAuelStandardForm(cubic : Nonflat);           
    else
        g, conicCoeffs := ConicFibrationForCpp(cubic);
    end if;
    discCoeffs := DiscriminantProjectionEquations(conicCoeffs);
        
    // If we have a cone, it is easier to count the points directly.
    if conicCoeffs[1..5] eq [0,0,0,0,0] then
        cc := PointCountsMagma(conicCoeffs[6] : Maxq:=Maxq);
        return [cc[i] * 2^(2*i) + 2^i + 1 : i in [1..Maxq]];
    end if;
    
    cppCounts := UncorrectedCppPointCounts(conicCoeffs, discCoeffs
                                           : ExecNum:=ExecNum, Minq:=Minq, Maxq:=Maxq,
                                             CompilerOptimization:=CompilerOptimization,
                                             CompileEachQ:=CompileEachQ,
                                             UseCache:=UseCache,
                                             PrintTimes:=PrintTimes);
    
    return CppCountCorrection(cppCounts, conicCoeffs);
end intrinsic;

intrinsic PointCounts(cubic, m) -> RngIntElt
{Compute the point counts on the given smooth cubic over F_(2^m). }
    return PointCounts(cubic : Minq:=m, Maxq:=m);
end intrinsic;

intrinsic PointCountsNoTransform(cubic : ExecNum:=false, Maxq:=11) -> SeqEnum
{Compute the pointcounts on the given cubic over Fq, where q = 2,4,...,2048.
In this version the cubic is assumed to be in Addington-Auel standard form.
Primarily, this function is used for testing.}

    conicCoeffs := ConicFibrationCoefficients(cubic);
    discCoeffs := DiscriminantProjectionEquations(conicCoeffs);

    // If we have a cone, it is easier to count the points directly.
    if conicCoeffs[1..5] eq [0,0,0,0,0] then
        cc := PointCountsMagma(conicCoeffs[6] : Maxq:=Maxq);
        return [cc[i] * 2^(2*i) + 2^i + 1 : i in [1..Maxq]];
    end if;
    
    cppCounts := UncorrectedCppPointCounts(conicCoeffs, discCoeffs
                                           : ExecNum:=ExecNum, Maxq:=Maxq);

    return CppCountCorrection(cppCounts, conicCoeffs);
end intrinsic;


intrinsic PointCountingSystemStrings(: ExecNum:=false,
                                       CompilerOptimization:=false,
                                       FixN:=false,
                                       UseCache:=true,
                                       Warn:=true) -> MonStgElt, MonStgElt, MonStgElt
{Format the compilation string, header file name, and the executable file name for calling g++.}

    compileString := "g++ -std=c++11";

    if not Warn then
        compileString *:= " -w";
    end if;

    // NOTE: On ARM chips, Magma emulates x86 anyways.
    processor := Read(POpen("arch", "r")); 
    if processor eq "arm\n" then
        compileString *:= " -DARM";
    else
        compileString *:= " -mpclmul";
    end if;
    
    if CompilerOptimization then
        compileString *:= " -O3";
    end if;

    if UseCache then
        compileString *:= " -DWITHCACHE";
    end if;
    
    // Option to ensure parallel code doesn't fight over executable file.
    if ExecNum cmpeq false then
	execFile := "a.out";
        headerFile := "coeffs.h";        
    else
	execFile := Sprintf("a.%o.out", ExecNum);
        headerFile := Sprintf("coeffs%o.h", ExecNum);
        compileString *:= Sprintf(" -DCOEFFSFILE='\"%o\"'", headerFile);
    end if;

    compileString *:= " -o " * execFile;

    // Finish compile string    
    if FixN cmpeq false then
        compileString *:= " tableio.cpp count.cpp";
    else
        compileString *:= Sprintf(" -DN=%o count_bigq.cpp", FixN);
    end if;
    
    return compileString, headerFile, execFile;
end intrinsic;


intrinsic CompileCppCounter(compileString : PrintTimes:=PrintTimes)
{}
    timeBeforeCompile := Cputime();
    retcode := System(compileString);
    if retcode ne 0 then
        error "Error at compile time step. Exit status: ", retcode;
    end if;

    if PrintTimes then
        print "Compile time: ", Cputime() - timeBeforeCompile;
    end if;

    return;
end intrinsic;


intrinsic UncorrectedCppPointCounts(conicCoeffs, discCoeffs
                                    : ExecNum:=false, Minq:= 1, Maxq:=11,
                                      CompilerOptimization:=false,
                                      CompileEachQ:=false,
                                      UseCache:=true,
                                      PrintTimes:=false) -> SeqEnum
{Call the C++ point counting engine. The output of this function is the correct point counts
*assuming* that the line defining the conic fibration lies in the smooth locus. Otherwise,
useful information is produced, but one has to correct the output.}

    // Change to C++ directory.
    entryDir := GetCurrentDirectory();
    ChangeDirectory(PATH_TO_LIB * "point_counting_cpp");
    
    if not CompileEachQ then
        // Create relevant strings.
        compileString, headerFile, execFile :=
            PointCountingSystemStrings( : ExecNum:=ExecNum,
                                          CompilerOptimization:=CompilerOptimization,
                                          UseCache:=UseCache);

    else
        // Just extract information for header files.
        _, headerFile, execFile := PointCountingSystemStrings( : ExecNum:=ExecNum);
    end if;

    // Prepare C++ code.
    ok_write := WriteHeaderFile(headerFile, conicCoeffs, discCoeffs);

    // Compile (if not fixing q)
    if not CompileEachQ then
        try
            CompileCppCounter(compileString : PrintTimes:=PrintTimes);
        catch e
            ChangeDirectory(entryDir);
            error e;
        end try;
    end if;
    
    // Execute.
    if PrintTimes then print "Run times:"; end if;
    cppOutputs := [];
    
    for m in [Minq .. Maxq] do

        if not CompileEachQ then
            processString := Sprintf("./%o %o", execFile, m);
        else
            warn := (m eq Minq);
            // Compile a version for q=2^m.
            compileString :=
            PointCountingSystemStrings( : ExecNum:=ExecNum,
                                          FixN:=m,
                                          CompilerOptimization:=true, // Force optimization.
                                          UseCache:=false, // No tables.
                                          Warn:=warn
                                      );

            try
                CompileCppCounter(compileString : PrintTimes:=PrintTimes);
            catch e
                ChangeDirectory(entryDir);
                error e;
            end try;

            processString := Sprintf("./%o", execFile);            
        end if;

        timeBeforeRun := Cputime();
        cppout := Read(POpen(processString, "r"));

        // Notify            
        if PrintTimes then
            print Sprintf("q=2^%o,", m), Cputime() - timeBeforeRun;
        end if;
        
        // Assign
        cppOutputs[m-Minq+1] := cppout;
    end for;
    
    try
	point_counts := [StringToInteger(out) : out in cppOutputs];
    catch e
        ChangeDirectory(entryDir);
	error "Error from C++. Received output: ", cppOutputs;
    end try;

    // Cleanup afterwards
    System(Sprintf("rm %o", execFile));

    if headerFile ne "coeffs.h" then
        System(Sprintf("rm %o", headerFile));
    end if;

    // Return to the current directory.
    ChangeDirectory(entryDir);
    return point_counts;
end intrinsic;

intrinsic PointCountsMagma(cubic : Maxq:=3) -> SeqEnum
{Compute the pointcounts on the given cubic over Fq, where q = 2,4,...,2048.
Primarily used for testing.
}
    R := Parent(cubic);
    P := Proj(R);
    X := Scheme(P, cubic);

    return [#RationalPoints(X, GF(2^i)) : i in [1..Maxq]];
end intrinsic;
