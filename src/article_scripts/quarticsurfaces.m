k := FiniteField(2);
G := GL(4, k);
R<[x]> := PolynomialRing(k, 4);
V, Bit := GModule(G, R, 4);
B := sub<V | Bit(x[1]^3*x[2] + x[2]^3*x[1])>;
U := V/B;
iota := Morphism(B, V);
pi := Morphism(V, U);
GU := ActionGroup(U);
phi := GModuleAction(U);
orbrepsU := [o[1] : o in Orbits(GU)];
stabs := [Stabilizer(GU, rep) : rep in orbrepsU];
orbreps := [];
time for i in [1..#orbrepsU] do

    v :=  orbrepsU[i] @@ pi;
    stab := stabs[i];
    Vstab, Bitstab := GModule(stab @@phi , R, 4);
    GVstab := ActionGroup(Vstab);
    coset := {v + iota(b) : b in B};

    // Set up the associative array.
    A := AssociativeArray();
    for u in coset do
        A[u] := false;
    end for;

    orbreps_i := [];
    time for u in coset do
        if A[u] eq true then continue; end if;
        orbit := Orbit(GVstab, u);
	orbitsize := #G/(#GVstab/#orbit);
        Append(~orbreps, <u, orbitsize>);
        // Update A to mark progress
        for g in orbit do
            A[g] := true;
        end for;
   end for;

end for;


//how many smooth quartic surfaces in P3? first remove zero orbit
orbreps := [o : o in orbreps | o[1] ne 0];
P := ProjectiveSpace(R);
&+[o[2] : o in orbreps | IsNonsingular(Scheme(P, o[1]@@Bit))];
//result of the above sum is 10590854400


//BONUS: cubic threefolds
N := 5; //number of variables (dim X = N-2, X in P^(N-1))
d := 3; // d=degree(X)

k := FiniteField(2);
G := GL(N, k);

R<[x]> := PolynomialRing(k, N);

V, bmap := GModule(GL(N,k), R, d);
P := ProjectiveSpace(R);
W, inclusion := sub<V | [bmap(x[i]^3) : i in [1..5]]>;
U := V/W;
pi := Morphism(V, U);
GU := ActionGroup(U);
phi := GModuleAction(U);

orbrepsU := [o[1] : o in Orbits(GU)];
stabs := [Stabilizer(GU, rep) : rep in orbrepsU];
orbreps := [];
for i in [1..#orbrepsU] do

    v :=  orbrepsU[i] @@ pi;
    stab := stabs[i];
    Vstab, Bitstab := GModule(stab @@phi , R, d);
    GVstab := ActionGroup(Vstab);
    coset := {v + inclusion(b) : b in W};

    // Set up the associative array.
    A := AssociativeArray();
    for u in coset do
        A[u] := false;
    end for;

    orbreps_i := [];
    time for u in coset do
        if A[u] eq true then continue; end if;
        orbit := Orbit(GVstab, u);
	orbitsize := #G/(#GVstab/#orbit);
        Append(~orbreps, <u, orbitsize>);
        // Update A to mark progress
        for g in orbit do
            A[g] := true;
        end for;
   end for;

end for;

//how many smooth cubic threefolds in P4? first remove zero orbit
orbreps := [o : o in orbreps | o[1] ne 0];
P := ProjectiveSpace(R);
&+[o[2] : o in orbreps | IsNonsingular(Scheme(P, o[1]@@bmap))];
//result of the above sum is 10239344640

