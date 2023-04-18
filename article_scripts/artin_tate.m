AttachSpec("../CubicLib/CubicLib.spec");
zetafunctions := ReadZetaFunctions();

// actualzetas[k] will be Q(t), the degree 22 factor in denominator of zeta function..
_<t> := PolynomialRing(Rationals());
actualzetas := AssociativeArray();
for k in Keys(zetafunctions) do
    f := Evaluate(zetafunctions[k], 4*t);
    if Coefficient(f,0) eq -1 then f := -f; end if;
    actualzetas[k] := f;
end for;

artintatevalues := AssociativeArray();
for k in Keys(zetafunctions) do
    f := actualzetas[k];
    // Indeed, I think Asher's email was not quite correct. Correcting the sign error
    // gives a value of 4 instead of 0.
    lvalue := Evaluate(actualzetas[k]/(1-4*t)^Valuation(f, 1-4*t), 1/4);

    // This is the special value of the zeta function ** divided by q^chi(X, OX, 2) **
    artintatevalues[k] := 1/16 * 8/9 * 1/lvalue;

    // The numerator is always 1, which is a good sign for sure.

    // Our guess is off by a factor of 9. Not sure what that's about.

    // H^{3,1} might collapse (by 1) in the supersingular case.
end for;

