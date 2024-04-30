AttachSpec("../CubicLib/CubicLib.spec");

differentiabilityDir := DatabaseDirectory();
smoothFile := differentiabilityDir * "smooth/smooth.csv";
singularFile := differentiabilityDir * "singular/singular.csv";

cubics := LoadCubicOrbitData(: Flat, Quick);

for i in [1..#cubics] do
    f := cubics[i];
    issmooth := IsSmooth(f);

    fname := IsSmooth(f) select smoothFile else singularFile;
    line := Sprintf("%o, %o", i, issmooth);
    Write(fname, line);
end for;
