// The following script draws the feasibility boxes reported in the paper.
AttachSpec("../CubicLib/CubicLib.spec");
TMO := "Too many orbits";
Y := true;
N := false;

function LastNonXEntry(n, allegedLastD, q)
    boo, tmo := IsFeasible(n, allegedLastD+1, q);
    assert tmo eq TMO;
    return IsFeasibleUnionFind(n, allegedLastD, q),
           IsFeasible(n, allegedLastD, q);
end function;

////////////////////
// q = 2
q := 2;

expectedLastEntries2 := [
    [Y, Y],
    [Y, Y],
    [N, Y],
    [Y, Y],
    [N, Y],
    [N, Y],
    [Y, Y],
    [Y, Y],
    [N, Y]
];

expectedLastD2 := [
    48, 8, 5, 3, 3, 3, 2, 2, 2    
];

for n in [1..9] do
    d := expectedLastD2[n];
    a, b := LastNonXEntry(n, d, q);
    assert [a,b] eq expectedLastEntries2[n];
end for;

// Cutoff for quadrics. I don't actually know when IsFeasible returns false, but
// it has gotten to the point where it takes too long to check.
if false then
    assert IsFeasible(21, 2, 2 : CheckOrbits:=false); // Takes forever...
end if;

////////////////////
// q = 3
q := 3;

expectedLastEntries3 := [
    [Y, Y],
    [N, Y],
    [N, Y],
    [N, Y],
    [N, N],
    [Y, Y],
    [N, N]    
];

expectedLastD3 := [
    31, 7, 4, 3, 3, 2, 2
];

for n in [1..7] do
    d := expectedLastD3[n];
    a, b := LastNonXEntry(n, d, q);
    assert [a,b] eq expectedLastEntries3[n];
end for;


////////////////////
// q = 5
q := 5;

expectedLastEntries5 := [
    [Y, Y],
    [N, Y],
    [Y, Y],
    [N, N],
    [Y, Y],
    [N, N]    
];

expectedLastD5 := [
    22, 6, 3, 3, 2, 2
];

for n in [1..6] do
    d := expectedLastD5[n];
    a, b := LastNonXEntry(n, d, q);
    assert [a,b] eq expectedLastEntries5[n];
end for;
