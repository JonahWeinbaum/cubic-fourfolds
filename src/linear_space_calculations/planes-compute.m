AttachSpec("../CubicLib/CubicLib.spec");

k := FiniteField(2);
lines := ReadLinesIndex();

// The 1395 echelon forms cutting out planes in P5.
planes := {@ EchelonForm(M) : M in Hom(VectorSpace(k, 3), VectorSpace(k, 6)) | Rank(M) eq 3 @};  

// A choice of four lines on each plane.
fourlines := [];
for i in  [1..#planes] do
    l := [ln : ln in lines | RowSpace(planes[i]) subset RowSpace(ln) ];
    fourlines[i] := {Index(lines, l[j]) : j in [1..4]};
end for;


function LoadLinesThroughCubics()
    // If the lines-compute.m script has been run, load the results.
    // i.e., load a list of lists `L` such that L[i][j] is the set of lines through the
    // (i,j)-th orbit representative (using the non-flattened indexing of orbit representatives).

    linesdir := DatabaseDirectory() * "linear_subspaces/lines_through_cubics/";
    echforms := ReadLinesIndex();

    result := [];
    for i in [1..85] do
        fname := Sprintf("lines-%o.data", i);
        fcontents := Read(linesdir * fname);

        // Parse the lists of indices.
        linesThroughLists := eval fcontents;

        // Obtain the corresponding echelon forms.
        subresult := [[echforms[ind] : ind in indices] : indices in linesThroughLists];

        // Append.
        Append(~result, subresult);
    end for;

    return result;
end function;

linesthrough := LoadLinesThoughCubics();

// The plane P goes through the i-th cubic if and only if four lines on that plane are contained
// in the i-th cubic. (The indexing is with respect to whatever indexing was used for
// linesthrough.)
planesthrough := [];
time for i in [1..#linesthrough] do
    planesthrough[i] := [];
    for j in [1..#planes] do
        if fourlines[j] subset linesthrough[i] then
            Append(~planesthrough[i], planes[j]);
        end if;
    end for;

    // Report.
    if (i mod 10000) eq 0 then print i; end if;
end for;
print "Incidence computation for planes finished.";


// Compute indices for more efficient storage.
planesthroughindexed := [];
planeindices := AssociativeArray();
for i in [1..#planes] do
    planeindices[planes[i]] := i;
end for;

for i in [1..#planesthrough] do
    planesthroughindexed[i] := [planeindices[p]: p in planesthrough[i]];
end for;

// Save the data.
SetColumns(0);
for i in [1..#planesthroughindexed] do
    datadir := DatabaseDirectory();
    subpath := "linear_subspaces/planes_through_cubics/";
    fname := Sprintf("planes-%o.data", i);
    Write(filename, planesthroughindexed[i]);
end for;



