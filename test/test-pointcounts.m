Attach("CubicLib.m");
load "computecharpoly.m";

//////////////////////////////////////////////
// Script start


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
function TestPointCounts(cubic)
    r<[x]> := PolynomialRing(FiniteField(2), 6);
    cubic := r ! cubic;
    Pcubic := ProjectiveSpace(Parent(cubic));
    X := Scheme(Pcubic, cubic);
    
    Maxq := 4;
    ptcts := PointCounts(cubic : Maxq := Maxq);
    ptcounts := [ptcts[i] : i in [1..Maxq]];

    // Expected values.
    verifiedPointCounts := [#Points(X, FiniteField(2^k)) : k in [1..Maxq]];
    
    assert verifiedPointCounts eq ptcounts;
    return true;
end function;


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

// Really long test -- 133 (s).
// time assert &and[TestPointCounts(f) : f in test_cubics];

print "Tests passed!";

