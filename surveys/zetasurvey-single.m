AttachSpec("../CubicLib/CubicLib.spec");

orbdata := LoadCubicOrbitData(: Flat:=true, Quick:=true); // 2 minute load.

for i in [1..20] do
    f := orbdata[i];
    try
	a := PointCounts(f);
	print a; // Want to parallelize this call.
	ok_write := WriteZetaData(i, IsSmooth(f), a);
    catch e
	print e;
    end try;
end for;
