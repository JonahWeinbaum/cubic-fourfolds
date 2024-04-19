
Attach("CubicLib.spec");

/*planes is the 1395 echelon forms cutting out planes in P5. fourlines is a choice of four lines on each cubic, while 

*/
k := FiniteField(2);
planes := {@ EchelonForm(M) : M in Hom(VectorSpace(k, 3), VectorSpace(k, 6)) | Rank(M) eq 3 @};  

fourlines := [];
lines := ReadLinesIndex();

for i in  [1..#planes] do
l := [ln : ln in lines | RowSpace(planes[i]) subset RowSpace(ln) ];
fourlines[i] := {Index(lines,l[1]), Index(lines,l[2]), Index(lines,l[3]), Index(lines,l[4])};
end for;

//Eventually replace with a .m script reading in serialized data.
load "readlinesinorbitsdata.m";

linesthrough := linesinorbitsdata;

planesthrough := [];
time for i in [1..#linesthrough] do
    planesthrough[i] := [];
    for j in [1..1395] do
        if fourlines[j] subset linesthrough[i] then
            planesthrough[i] :=  planesthrough[i] cat [planes[j]];
        end if;
    end for;
    if (i mod 10000) eq 0 then print i;
    end if;
end for;
print "data computation finished";

planesthroughindexed := [];
planeindices := AssociativeArray();
for i in [1..#planes] do
planeindices[planes[i]] := i;
end for;

for i in [1..#planesthrough] do
planesthroughindexed[i] := [planeindices[p]: p in planesthrough[i]];
end for;

filename := "../..//database/linear_subspaces/planes_through_cubics/nonserialized-planes-data.csv";
SetColumns(0);
time for i in [1..#planesthroughindexed] do
Write(filename, planesthroughindexed[i]);
end for;
print "data file writing finished";


