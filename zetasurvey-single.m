Attach("CubicLib.m");

// Load Point counting function.
load "computecharpoly.m";
orbdata := LoadCubicOrbitData(: Flat:=true); // 2 minute load.

for i in [1..20] do
    b := orbdata[i]; // For example;

    try
	f := BitListToCubic(b);
	a := PointCounts(f);
	print a; // Want to parallelize this call.
	WriteZetaData(i, IsSmooth(f), a);
    catch e
	print e;
    end try;
end for;



