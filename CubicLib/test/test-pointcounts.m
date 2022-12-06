AttachSpec("CubicLib.spec");
//load "computecharpoly.m";

//////////////////////////////////////////////
// Script start

function CheckSpecialPointIsSingular(conicCoeffs, v)
    A, B, C, D, E, F := Explode(conicCoeffs);
    discP3 := A*E^2 + B^2*F + C*D^2 - B*D*E;
    X := Scheme(Proj(Parent(discP3)), discP3);
    assert Dimension(X) eq 2;

    seqv := Eltseq(v);
    assert Evaluate(B, seqv) eq 0;
    assert Evaluate(D, seqv) eq 0;
    assert Evaluate(E, seqv) eq 0;
    
    assert IsSingular(X ! Eltseq(v));
    return true;
end function;


// Checks that the conversion to Addington-Auel standard form works.
function TestAAStdForm(cubic)
    r<[x]> := PolynomialRing(FiniteField(2), 6);
    cubic := r ! cubic;
    g, coeffs := AddingtonAuelStandardForm(cubic);
    P := ProjectiveSpace(Parent(g));

    assert IsSubscheme(Scheme(P, [x[1], x[2], x[3], x[4]]), Scheme(P, g));
    assert CheckSpecialPointIsSingular(coeffs, [0,0,0,1]);
    
    return true;
end function;


// Checks if C++ point counting matches magma's point counts.
function TestPointCounts(cubic : ExecNum:=false, PrintTimes:=false)
    r<[x]> := PolynomialRing(FiniteField(2), 6);
    cubic := r ! cubic;
    Pcubic := ProjectiveSpace(Parent(cubic));
    X := Scheme(Pcubic, cubic);
    
    Maxq := 4;
    ptcts := PointCounts(cubic : Maxq := Maxq, ExecNum:=ExecNum, PrintTimes:=PrintTimes);
    ptcounts := [ptcts[i] : i in [1..Maxq]];

    // Expected values.
    verifiedPointCounts := [#Points(X, FiniteField(2^k)) : k in [1..Maxq]];
    
    assert verifiedPointCounts eq ptcounts;
    return true;
end function;


////////////////////////////
// Tests Begin

// Ambients
P5<[x]> := ProjectiveSpace(GF(2), 5);

load "test/test-cubics.m"; // Gives us a list of cubics.

print "Beginning tests...";

print "Addington-Auel standard form...";
time for f in test_cubics do
    assert TestAAStdForm(f);
end for;

print "Point counts...";
time for f in test_cubics[1..5] do
    assert TestPointCounts(f);
end for;

print "Point counts (with options)...";
time for f in test_cubics[1..5] do
    assert TestPointCounts(f : ExecNum:=1);
end for;

print "Benchmarks...";
time for f in test_cubics[6..6] do
         dummy := PointCounts(f : PrintTimes);
end for;


// Really long test -- 133 (s).
// time assert &and[TestPointCounts(f) : f in test_cubics];

print "Tests passed!";

