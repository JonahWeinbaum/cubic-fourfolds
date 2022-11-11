SetColumns(0);

//////////////////////////////////////////////
// Script start
AttachSpec("CubicLib.spec");
//load "data-processing/read-data-lines-index.m";
//load "data-serialization/dataprocessing-lines.m";



if not assigned COMPUTE_CHARPOLY_ALREADY_LOADED then
    // This block of code is used to set up the LinesThrough
    // function. Once it loads properly, 
        
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

        monoevaluated[form][1] := BitwiseOr(ShiftLeft(monoevaluated[form][1], 1), Integers()!eval1);
        monoevaluated[form][2] := BitwiseOr(ShiftLeft(monoevaluated[form][2], 1), Integers()!eval2);
        monoevaluated[form][3] := BitwiseOr(ShiftLeft(monoevaluated[form][3], 1), Integers()!eval3bas1);
        monoevaluated[form][4] := BitwiseOr(ShiftLeft(monoevaluated[form][4], 1), Integers()!eval3bas2);
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

// Given a cubic polynomial f, return the set of lines through f.
function LinesThrough(f)

    // Type conversion back to bitstring.
    //b := Reverse([Integers() ! Bit(f)[i] : i in [1..56]]);
    //f := Seqint(b, 2);

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
end function;

COMPUTE_CHARPOLY_ALREADY_LOADED := true;
end if;


function AddingtonAuelStandardForm(cubic : Nonflat:=false)
    // Transforms a cubic into a form as outlined in Addington-Auel Section 3.
    
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
end function;


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

    Write(fname, string : Overwrite);
    return 0;
end function;
    
function PointCounts(cubic : ExecNum:=false, Maxq:=11, Nonflat:=false)
    // This function is mostly based off of 
    // Section 3 of Addington-Auel. See loc. cit. for details. 

    k := GF(2);
    V6 := VectorSpace(k, 6);
    R<y0, y1, y2, y3> := PolynomialRing(k, 4);
    RR<y4, y5> := PolynomialRing(R, 4);

    // Now compute coefficients for the conic bundle fibration
    // that will be fed into C++ point counting code.
    if Nonflat then
     g, conicCoeffs := AddingtonAuelStandardForm(cubic : Nonflat:=true);
     else   g, conicCoeffs := AddingtonAuelStandardForm(cubic);
    end if;
    A, B, C, D, E, F := Explode(conicCoeffs);
    a, b, c, d := DiscriminantProjectionEquations(conicCoeffs);
    discCoeffs := [a,b,c,d];

    
    // Testing code.
    /* X := Scheme(Proj(Parent(cubic)), cubic); */
    /* print [#Points(X, GF(2^j)) : j in [1..4]]; */

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
    ChangeDirectory("point_counting_cpp");
    
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
end function;


/*
// Zeta function stuff.
function zetaFun()
    // Stuff for zeta function.
    tr := [(StringToInteger(Read(POpen("./a.out " cat Sprint(m), "r"))) - 1 - 2^m - 4^m - 8^m - 16^m)/4^m : m in [1..11]];
    
    c := [];
    c[1] := -tr[1];
    for k in [2..11] do
	c[k] := (tr[k] + &+[c[i]*tr[k-i] : i in [1..k-1]] )/(-k);
    end for;

    p := t^22 + &+[c[i]*t^(22-i) : i in [1..11] ];
    g := QQt ! (t^22 * Evaluate(p, 1/t));

    // TODO: fix this.
    if c[11] ne 0 then
	ret := p + Modexp(g, 1, t^11);
    else
	ret := p - Modexp(g, 1, t^11);
    end if;

    return;
end function;
*/
