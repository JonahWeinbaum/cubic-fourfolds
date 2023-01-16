// Setup jobs for the Discovery cluster.
AttachSpec("CubicLib/CubicLib.spec");

function OKTransform(f)
    try
        PrepareClusterCppComputation(f);
        return true;
    catch e
        return false;
    end try;
end function;


// Test cubics.
P5<[x]> := ProjectiveSpace(GF(2), 5);
load "CubicLib/test/test-cubics.m";
cubics := [f : f in test_cubics | OKTransform(f)];


// 1. Choose the problematic cubics.

counts := ReadCSV("zeta_functions/point_counts.csv");
cubics := LoadCubicOrbitData(: Flat);

// Filter based on how many zeros are there.




// 2. Prepare jobs.
for i in [1..#cubics] do
    f := cubics[i];
    PrepareClusterCppComputation(f : ExecNum:=i);
end for;

// 3.(??) Copy to Discovery?
