AttachSpec("../CubicLib/CubicLib.spec");

print "Running in Quick mode. We *strongly* recommend running the parallel version." *
"If you really want to run the serialized version, change Quick to false in the script file.";

Quick := true; // Change this to false if you don't mind waiting a few days.

orbdata := LoadCubicOrbitData(: Flat, Quick:=Quick);

for i in [1..#orbdata] do
    f := orbdata[i];
    try
	a := PointCounts(f);
	print a; // Want to parallelize this call.
	ok_write := WriteZetaData(i, IsSmooth(f), a);
    catch e
	print e;
    end try;
end for;
