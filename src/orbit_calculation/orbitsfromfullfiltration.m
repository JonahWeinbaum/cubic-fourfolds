/// We'll use the filtration V > Wp > W to compute a complete set of orbit representatives for G acting on V

load "dataprocessing.m";

k := FiniteField(2);
G := GL(6, k);
P5<[x]> := ProjectiveSpace(k, 5);
R := CoordinateRing(P5);


V, Bit := GModule(G, R, 3);
W      := sub<V | Bit(x[1]^3) > ;
U, piU := quo<V | W>;

phi := Morphism(W, V);

forms, formmap := GModule(G, R, 1);

GV := ActionGroup(V);
GW := ActionGroup(W);
GU := ActionGroup(U);

//Hom from G to GL(U)
umap := GModuleAction(U);

//space of quadratic*linear, which has dimension 36
Wp := sub<V | Bit(x[1]^2*x[2])>;


Up, phip := V/Wp;
upmap := GModuleAction(Up);
GUp := ActionGroup(Up);

// Construct orbit representatives for the first orbit space.
GUporbits := Orbits(GUp);
GUporbreps := [GUporbits[i][1] : i in [1..#GUporbits]]; 

// Also determine the stabilizers and lift back to the original group.
stabsUp := [Stabilizer(GUp, a) : a in GUporbreps];
stabsUpinG := [H @@ upmap : H in stabsUp];

///////////////////////
// Find V/W orbit reps

orbvwreps := [];
stabsvw := [];
j:= 0;

for i in [1..#stabsUpinG] do

    rep := GUporbreps[i];
    stab := stabsUpinG[i];

    // Compute the restriction module.
    Vstab, Bitstab := GModule(stab, R, 3);   
    Wstab := sub<Vstab | [Vstab ! Vector(phi(w))  : w in Basis(W)]>;
    Ustab, pistab := Vstab/Wstab;
    ustabmap := GModuleAction(Ustab);
    GUstab := ActionGroup(Ustab);

    coset := {piU(rep @@phip) + w : w in piU(Wp)};

    // Set up the associative array.
    A := AssociativeArray();
    for v in coset do
        A[v] := false;
    end for;

    time for v in coset do
             if A[v] eq true then continue; end if;
             orbit := Orbit(GUstab, v);
             Append(~orbvwreps, v);
	     j := j+1;
             // Update A to mark progress
             for g in orbit do
                 A[g] := true;
             end for;

             print j;
             // Record the stabilizer of the representative for use higher up in the filtration.
             stabsvw[j] := Stabilizer(GUstab, orbvwreps[j]) @@ustabmap;
         end for;
    
end for;


print "Precomputation 1: complete.";
print "Starting orbit computation...";


for i in [1..85] do
	
    fname := Sprintf("orbitreps-%o.data", i);
    file := Open(fname, "w");

    // Construct the restriction G-module
    v := orbvwreps[i];
    stab := stabsvw[i];
    Vstab, Bitstab := GModule(stab, R, 3);
    GVstab := ActionGroup(Vstab);
    coset := {(v @@ piU) + phi(w) : w in W};

    // Set up the associative array.
    A := AssociativeArray();
    for u in coset do
        A[u] := false;
    end for;

    orbreps_i := [];
    time for u in coset do
        if A[u] eq true then continue; end if;
        orbit := Orbit(GVstab, u);
        Append(~orbreps_i, u);
	WriteBytes(file, serialize(u));

        // Update A to mark progress
        for g in orbit do
            A[g] := true;
        end for;

    end for;

    print i, #orbreps_i;
    print "Computation ", i, "/85 complete.";
end for;

