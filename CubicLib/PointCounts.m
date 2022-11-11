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
    yvec := Vector(R, [y0,y1,y2,y3]);
    newcoeffs := [Evaluate(co, Eltseq(yvec * ChangeRing(invM4, R))) : co in coeffs];

    return g, newcoeffs;
end intrinsic;


function DiscriminantProjectionEquations(conicCoeffs)
    A, B, C, D, E, F := Explode(conicCoeffs);

    rrr<y0, y1, y2> := PolynomialRing(k, 3);
    RRR<y3> := PolynomialRing(rrr, 1);
    disc := Evaluate(C*D^2 + A*E^2 + F*B^2 + B*D*E, [RRR!y0, RRR!y1, RRR!y2, RRR!y3]);

    a := MonomialCoefficient(disc, y3^3);
    b := MonomialCoefficient(disc, y3^2);
    c := MonomialCoefficient(disc, y3);
    d := MonomialCoefficient(disc, 1);    
    
    return a, b, c, d;
end function;
    

function WriteHeaderFile(fname, conicCoeffs, discCoeffs)

    A, B, C, D, E, F := Explode(conicCoeffs);
    a, b, c, d := Explode(discCoeffs);
    
    string :=  "#define abcd \\" cat "\n" cat 
               "    a = " cat CppInputString(a) cat ", \\" cat "\n" cat
	       "    b = " cat CppInputString(b) cat ", \\" cat "\n" cat
	       "    c = " cat CppInputString(c) cat ", \\" cat "\n" cat
	       "    d = " cat CppInputString(d) cat "\n" cat "\n" cat
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

intrinsic PointCounts(cubic : ExecNum:=false, Maxq:=11, Nonflat:=false) -> SeqEnum
{Compute the pointcounts on the given cubic over Fq, where q = 2,4,...,2048.}

    // This function is mostly based off of 
    // Section 3 of Addington-Auel. See loc. cit. for details.
    //
    // The actual driver of this function is the C++ code in point_counting_cpp.

    k := GF(2);
    V6 := VectorSpace(k, 6);
    R<y0, y1, y2, y3> := PolynomialRing(k, 4);
    RR<y4, y5> := PolynomialRing(R, 4);

    // Now compute coefficients for the conic bundle fibration
    // that will be fed into C++ point counting code.
    if Nonflat then
        g, conicCoeffs := AddingtonAuelStandardForm(cubic : Nonflat:=true);
    else
        g, conicCoeffs := AddingtonAuelStandardForm(cubic);
    end if;
    
    A, B, C, D, E, F := Explode(conicCoeffs);
    a, b, c, d := DiscriminantProjectionEquations(conicCoeffs);
    discCoeffs := [a,b,c,d];

    // Option to ensure parallel code doesn't fight over executable file.
    if ExecNum cmpeq false then
	execFile := "a.out";
        headerFile := "coeffs.h";
        compileString := Sprintf("g++ tableio.cpp count.cpp -o %o", execFile);
    else
	execFile := Sprintf("a.%o.out", ExecNum);
        headerFile := Sprintf("coeffs%o.h", ExecNum);
        compileString := Sprintf("g++ -DWITHCACHE -DCOEFFSFILE='\"%o\"' tableio.cpp count.cpp -o %o", headerFile, execFile);
    end if;

    // Change to C++ directory
    entryDir := GetCurrentDirectory();
    ChangeDirectory(PATH_TO_LIB * "point_counting_cpp");
    
    // Prepare C++ code.
    ok_write := WriteHeaderFile(headerFile, conicCoeffs, discCoeffs);

    // Compile.
    retcode := System(compileString);
    error if retcode ne 0, "Error at compile time step. Exit status: ", retcode;
    
    // Execute.
    cppOutputs := [Read(POpen("./" * execFile * " " cat Sprint(m), "r")) : m in [1..Maxq]];

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

