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


function LoadLinesThroughCubics( : ReturnIndices:=false)
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
        if ReturnIndices then
            subresult := [Set(indices) : indices in linesThroughLists];
        else
            subresult := [{echforms[ind] : ind in indices} : indices in linesThroughLists];
        end if;
        
        // Append.
        Append(~result, subresult);
    end for;

    return result;
end function;

linesthrough := LoadLinesThroughCubics(:ReturnIndices);
assert #linesthrough eq 85;

// The plane P goes through the (i,j)-th cubic if and only if four lines on that plane are
// contained in the (i,j)-th cubic. (The indexing is with respect to the usual (i,j)-indexing
// for orbit representatives obtained from the filtration, or equivalently,
// the (i,j)-th cubic is LoadCubicOrbitData(: Flat:=false)[i][j];

// File info.
SetColumns(0);
datadir := DatabaseDirectory();
subpath := "linear_subspaces/planes_through_cubics/";

for i in [1..#linesthrough] do

    subplanesthrough := [];
    for j in [1..#linesthrough[i]] do
        subplanesthrough[j] := [];

        // Check if each plane is contained in the (i,j)-th cubic via the lines.
        // If so, save the index.
        for ll in [1..#planes] do
            if fourlines[ll] subset linesthrough[i][j] then
                Append(~subplanesthrough[j], ll);
            end if;
        end for;
    end for;

    // Write to file.
    fname := Sprintf("planes-%o.data", i);
    Write(datadir * subpath * fname, subplanesthrough);
    
    // Report.
    print Sprintf("Completed %o / %o tasks.", i, #linesthrough);
end for;
print "Incidence computation for planes finished.";



