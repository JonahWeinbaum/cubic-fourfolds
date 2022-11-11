AttachSpec("../CubicLib/CubicLib.spec");

// Load Point counting function.
orbdata := LoadCubicOrbitData(: Flat, Quick); // 2 minute load.

for i in [1..20] do
    b := orbdata[i]; // For example;

    try
	f := BitListToCubic(b);
	a := CubicOrbitSize(f);
	print a; // Want to parallelize this call.
	ok_write := WriteOrbitSizeData(i, a);
    catch e
	print e;
    end try;
end for;



