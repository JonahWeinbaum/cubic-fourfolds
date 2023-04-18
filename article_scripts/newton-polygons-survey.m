AttachSpec("../CubicLib/CubicLib.spec");
zetafunctions := ReadZetaFunctions();

//actualzetas[k] will be Q(t), the degree 22 factor in denominator of zeta function..
_<t> := PolynomialRing(Rationals());

actualzetas := AssociativeArray();
for k in Keys(zetafunctions) do
    f := Evaluate(zetafunctions[k], 4*t);
    if Coefficient(f,0) eq -1 then f := -f; end if;
    actualzetas[k] := f;
end for;

nps := AssociativeArray();

for k in Keys(actualzetas) do
    nps[k] := LowerVertices(NewtonPolygon(actualzetas[k], 2));
end for;


//Note: the hodge polygon has vertices [ <0, 0>, <1, 1>, <21, 41>, <22, 44> ];
hodgepolygon := [ <0, 0>, <1, 1>, <21, 41>, <22, 44>];

