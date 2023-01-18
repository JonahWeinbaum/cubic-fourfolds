// Setup jobs for the Discovery cluster.
AttachSpec("CubicLib/CubicLib.spec");
load "CubicLib/Fano.m";

function OKTransform(f)
    try
        PrepareClusterCppComputation(f);
        return true;
    catch e
        return false;
    end try;
end function;

function BadnessOfCounts(counts)
    return Badness(HalfWeil(counts));
end function;


// Test cubics.
// P5<[x]> := ProjectiveSpace(GF(2), 5);
// load "CubicLib/test/test-cubics.m";
// cubics := [f : f in test_cubics | OKTransform(f)];


// 1. Choose the problematic cubics.

time counts := ReadCSV("zeta_functions/point_counts.csv");
time cubics := LoadCubicOrbitData(: Flat);


// Set up dictionaries
time cubicsDic := InverseAssociativeArray(SequenceAsDictionary(cubics) : UniqueKeys);
time cubcounts := Compose(cubicsDic, counts : IgnoreMissingKeys);

// Filter by badness
badness_filter := func<C | BadnessOfCounts(C) eq 10>;
time cubcounts_bad10 := CacheFilter(badness_filter, cubcounts);

// Filter by whether the Fano scheme is empty.
relevant := AssociativeArray();
for k in Keys(cubcounts_bad10) do
    if Dimension(FanoSchemePatch(k)) eq -1 then
        relevant[k] := cubcounts_bad10[k];
    end if;
end for;


// 2. Prepare jobs.
for f in Keys(relevant) do
    try
        PrepareClusterCppComputation(f : ExecNum:=Index(cubics, f));
    catch e
        print "Error at:", Index(cubics, f);
        print e;
        // do_nothing := true;
    end try;
end for;

// TODO: XXX: Apparently the change-of-coordinate rank issues occasionally pop up.
//            meaning I should actually find some kind of extra algebraic cycle?

// 3.(??) Copy to Discovery?
//   NOTE: Requires a VPN connection to Dartmouth's network.
