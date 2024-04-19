AttachSpec("CubicLib.spec");

// Ambients
P5<[x]> := ProjectiveSpace(GF(2), 5);
load "test/test-cubics.m"; // Gives us a list of cubics.
print "Beginning tests...";


// We compare the big C++ point counts to the actual values calculated from the
// zeta function. (In cases where we know the sign of the functional equation.)


// This test takes a long time to run.
ambiguousList    := [4, 5, 9, 25, 28, 33, 39, 54, 67, 68, 70, 73, 80, 94, 95];
badTransformList := [1, 7, 21, 24, 30, 49, 64, 66, 79, 82, 85, 88, 91];
bugged_RankDegen := [29, 43, 58, 61, 62, 63, 90, 96, 97];


// NOTE: Non-surjectivity of the projection *might* be a problem, but I don't think
// it can happen for smooth cubics.

B := 5; // Pointcount bound.
for i in [1..100] do
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
    new_counts := PointCounts(f : Minq := 1, Maxq := B, CompileEachQ:=true);

    if not extended_counts[1..B] eq new_counts then
        ec := extended_counts[1..B];
        nc := new_counts[1..B];
        error "Mismatched point counts:", [ec[i] - nc[i] : i in [1..B]];
    end if;
    
end for;


// Test a specific value.
function testSpecific(i)
    f := test_cubics[i];
    counts := PointCounts(f);
    chi := Charpoly(counts);
    extended_counts := CubicWeilPolynomialToPointCounts(chi);

    // First test the conversion function.
    assert extended_counts[1..11] eq counts;

    // Next, test the big point counts.
    new_counts := PointCounts(f : Minq := 1, Maxq := B, CompileEachQ:=true);

    if not extended_counts[1..B] eq new_counts then
        ec := extended_counts[1..B];
        nc := new_counts[1..B];
        error "Mismatched point counts:", [ec[i] - nc[i] : i in [1..B]];
    end if;
    return 0;
end function;
