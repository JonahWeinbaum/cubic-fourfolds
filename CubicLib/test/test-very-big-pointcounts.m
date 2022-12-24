AttachSpec("CubicLib.spec");

// Ambients
P5<[x]> := ProjectiveSpace(GF(2), 5);
load "test/test-cubics.m"; // Gives us a list of cubics.
print "Beginning tests...";


// We compare the big C++ point counts to the actual values calculated from the
// zeta function. (In cases where we know the sign of the functional equation.)

// WARNING: This test takes *ages* to run.
ambiguousList    := [4, 5, 9, 25, 28, 33, 39, 54, 67, 68, 70, 73, 80, 94, 95];
badTransformList := [1, 7, 21, 24, 30, 49, 64, 66, 79, 82, 85, 88, 91];
bugged_RankDegen := [29, 43, 58, 61, 62, 63, 90, 96, 97];

B := 18;
for i in [20..28] do
    if i in ambiguousList    then continue; end if;
    if i in badTransformList then continue; end if;
    // TODO: Might be nice to handle these cases.
    if i in bugged_RankDegen then continue; end if;

    print i;
    f := test_cubics[i];
    counts := PointCounts(f);
    chi := Charpoly(counts);
    extended_counts := CubicWeilPolynomialToPointCounts(chi);
    
    // First test the conversion function.
    assert extended_counts[1..11] eq counts;

    // Next, test the big point counts.
    new_counts := PointCounts(f : Minq := B, Maxq := B, CompileEachQ:=true);
    assert extended_counts[B] eq new_counts[1];
end for;

// ~  1 minute  for 2^15
// ~  4 minutes for 2^16
// ~ 16 minutes for 2^17
// ~  1 hour    for 2^18
