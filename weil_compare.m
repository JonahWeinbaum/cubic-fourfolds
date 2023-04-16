AttachSpec("CubicLib/CubicLib.spec");

////////////////////////////////////////////////////////////////////////////////
// Task 1 : Count Kedlaya-Sutherland objects and check out numbers match theirs.

fileLines := ReadLines(DatabaseDirectory() * "zeta_functions/kedlaya_sutherland_list.csv");
transPolys := [Polynomial(eval(line)) : line in fileLines];

// To count properly, we need to count the cyclotomic polynomials of each degree
// and count partitions.

cyclo := [f : f in [CyclotomicPolynomial(d) : d in [1..100]] | Degree(f) le 21];

_<T> := PowerSeriesRing(Rationals(), 22);
comboSeries := &*[1/(1-T^Degree(f)) : f in cyclo];

totalNum := &+[Coefficient(comboSeries, 21-Degree(f)) : f in transPolys];
assert totalNum eq 2971182;


////////////////////////////////////////////////////////////////////////////////
// Task 2 : Count the number of zeta functions from our database.
ourWeilPolys := eval Read(DatabaseDirectory()*"zeta_functions/zeta_coeffs_NEWandCORRECT.csv");

weilPolySet := {Polynomial(tup[2]) : tup in ourWeilPolys};
assert #weilPolySet eq 86472;


////////////////////////////////////////////////////////////////////////////////
// Task 3 : Count the number of zeta functions from our database.

function KSNormalization(f)
    return 2 * f / Evaluate(f, 0);
end function;

function Theorem2Test(f)
    R<t> := Parent(f);

    // Make sure that f is normalized to match the
    // conditions of Kedlaya-Sutherland
    L := KSNormalization(f);
    L1 := L/(1-t)^Valuation(L, 1-t);

    return IsSquare(Evaluate(L1, -1));
end function;


function Computation3cTest(f)
    q := 2;
    R<t> := Parent(f);
    L := KSNormalization(f/(1-t));

    K<T> := PowerSeriesRing(Rationals(), 10);
    den := (1 - T) * (1 - q*T) * (1 - q^2 * T) * q^(-1) * Evaluate(L, q*T);

    
    try
        series := -Log(den);
    catch e
        print den;    
        error e;
    end try;

    // Check for the conditions.
    assert Coefficient(series, 0) eq 0;
    coeffs := [Coefficient(series, i) : i in [1..9]];

    allPositive := &and [c ge 0 : c in coeffs];
    if not allPositive then return false; end if;

    ksCondition := true;
    for tup in [<2,1>, <3, 1>, <4, 2>] do
        nm := tup[1]; n := tup[2];
        ksCondition and:= (nm * coeffs[nm] ge n * coeffs[n]);
    end for;

    return ksCondition;
end function;


function BothTests(f)
    return Theorem2Test(f) and Computation3cTest(f);
end function;

// We need at least one non-Lefschetz algebraic cycle over F2 for any of this to make sense.
R<t> := Parent(Random(weilPolySet));
specialPolySet := {f : f in weilPolySet | Valuation(f, 1-t) gt 0};

inKSlist := {f : f in specialPolySet | BothTests(f)};
notinKSlist := {f : f in specialPolySet | not BothTests(f)};

