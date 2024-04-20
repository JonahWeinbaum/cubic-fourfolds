AttachSpec("../CubicLib/CubicLib.spec");

print "Running in Quick mode. We *strongly* recommend running the parallel version." *
"If you really want to run the serialized version, change Quick to false in the script file.";

Quick := true; // Change this to false if you don't mind waiting a few days.

// Load Point counting function.
orbdata := LoadCubicOrbitData(: Flat, Quick:=Quick);

for i in [1..#orbdata] do
    f := orbdata[i];

    try
	a := CubicOrbitSize(f);
	ok_write := WriteOrbitSizeData(i, a);
    catch e
	print e;
    end try;
end for;



