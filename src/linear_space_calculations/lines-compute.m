////////////////////////////////////////////////////////////////////////////////
//
// lines-compute.m
//
////////////////////////////////////////////////////////////////////////////////
//
// This file finds the list of lines through every orbit representative
// within our database. The result is a CSV file.

AttachSpec("../CubicLib/CubicLib.spec");

// Read the collection of lines in the ambient projective space.
lines := ReadLinesIndex();

k := FiniteField(2);
F4 := FiniteField(4);
basF4 := Basis(F4);
k6 := RSpace(k, 6);
k4 := RSpace(k, 4);

G_4 := GL(6, F4);
P5_4<[x]> := ProjectiveSpace(F4, 5);
R_4 := CoordinateRing(P5_4);

V_4, Bit_4 := GModule(G_4, R_4, 3);

ptsonlines := AssociativeArray();
monoevaluated := AssociativeArray();

// Loop over all lines.
for form in lines do
    nullsp := NullspaceOfTranspose(form);
    bas := Basis(nullsp);
    nullsp := ExtendField(nullsp, F4);
    pt1 := bas[1];
    pt2 := bas[2];
    pt3 := (nullsp!pt1)*1 + (nullsp!pt2)*F4.1;

    ptsonlines[form] := <pt1, pt2, pt3>;

    // Allocate and populate the 4-tuple of evaluations.
    monoevaluated[form] := <0,0,0,0>;
    for mono in Basis(V_4) do 
        eval1 := Integers() ! Evaluate(mono @@ Bit_4, Eltseq(pt1));
        eval2 := Integers() ! Evaluate(mono @@ Bit_4, Eltseq(pt2));
        
        eval3 := Evaluate(mono @@ Bit_4, Eltseq(pt3));
        eval3bas1 := Integers() ! Eltseq(eval3)[1];
        eval3bas2 := Integers() ! Eltseq(eval3)[2];

        monoevaluated[form][1] := BitwiseOr(ShiftLeft(monoevaluated[form][1], 1), eval1);
        monoevaluated[form][2] := BitwiseOr(ShiftLeft(monoevaluated[form][2], 1), eval2);
        monoevaluated[form][3] := BitwiseOr(ShiftLeft(monoevaluated[form][3], 1), eval3bas1);
        monoevaluated[form][4] := BitwiseOr(ShiftLeft(monoevaluated[form][4], 1), eval3bas2);
    end for;
end for;

// Compute the incidence relations of cubics and lines.
orbdata := LoadCubicOrbitData(: Flat:=false);

k := FiniteField(2);
G := GL(6, k);
P5<[x]> := ProjectiveSpace(k, 5);
R := CoordinateRing(P5);
V, Bit := GModule(G, R, 3);

// Compute lines through all cubics in database. Store the indices for more efficient loading.
for i in [1..85] do
    o := orbdata[i];
    linedata := [LinesThrough(f : ReturnIndices) : f in o];
    print i;

    // Write to file.
    datadir := DatabaseDirectory();
    subpath := "linear_subspaces/lines_through_cubics/";
    fname := Sprintf("lines-%o.data", i);
    PrintFile(datadir * subpath * fname, linedata);
end for;
