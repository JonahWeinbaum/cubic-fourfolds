//I. Precomputation
AttachSpec("CubicLib.spec");

/*planes is the 1395 echelon forms cutting out planes in P5. fourlines is a choice of four lines on each cubic
*/

k := FiniteField(2);
planes := {@ EchelonForm(M) : M in Hom(VectorSpace(k, 3), VectorSpace(k, 6)) | Rank(M) eq 3 @};

lines := ReadLinesIndex();

//we are distinguishing planes by choosing 4 lines on each one
fourlines := [];
for i in [1..#planes] do
l := [ln : ln in lines | RowSpace(planes[i]) subset RowSpace(ln) ];
fourlines[i] := {l[1], l[2], l[3], l[4]};
end for;


//II. PlanesThrough function

function PlanesThrough(cubic)
  linesthrough := LinesThrough(cubic);
   planesthrough := [];
   for j in [1..1395] do
if fourlines[j] subset linesthrough then
planesthrough :=planesthrough cat [planes[j]];
end if;
end for;
return planesthrough;
end function;
