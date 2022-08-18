///////////////////////////////////
// cubics.m

k := GF(2);
P5<[x]> := ProjectiveSpace(k, 5);
R := CoordinateRing(P5);

// Compute the norms  of linear forms
function Norm(ll)
    d := Degree(BaseRing(Parent(ll)));
    coeffs, mons := CoefficientsAndMonomials(ll);

    Cs := [[c^(2^i) : c in coeffs] : i in [0..d-1]];
    Ls := [Polynomial(C, mons) : C in Cs];
    
    return &*Ls;
end function;

// Count Orbits using Burnside.
function CountOrbits(Gmod)
    rho := Representation(Gmod);
    G := Group(Gmod);
    return &+[#Eigenspace(rho(c[3]),1) * c[2] : c in ConjugacyClasses(G) ]/#G;
end function;

G := GL(6, 2);

V, toV := GModule(G, R, 3);
O1, phi :=  GModule(G, R, 1);
iphi := Inverse(phi);


//////////////////
// Squares times linears

square_times_linear := {iphi(l1) * iphi(l2)^2 : l1 in O1, l2 in O1};
W1 := sub<V | [toV(R ! g) : g in square_times_linear]>;

// Burnside computation.
Vmod := quo<V | W1>;
num_orbits := CountOrbits(Vmod);

/////////////////////////
// k-forms that become Warring representable over extensions.

K := GF(2^3);
GK := ChangeRing(G, K);
P5K := ProjectiveSpace(K, 5);
RK := CoordinateRing(P5K);

O1K, phiK := GModule(GK, RK, 1);
iphiK := Inverse(phiK);
VK, toVK := GModule(GK, RK, 3);

// Warring representable over K
cubes := {iphiK(ll)^3 : ll in O1K};
WK := sub<VK | [toVK(g) : g in cubes]>;

ResVK := VectorSpace(K, Degree(K) * Dimension(VK));
RatSubspace := sub<ResVK | [ResVK.i : i in [1..Dimension(ResVK) by Degree(K)]]>;

function WeilRestrict(g)
    h := VK ! g;
    hvec := &cat [Eltseq(a) : a in Eltseq(h)];

    return ResVK ! hvec;
end function;

ResWK := sub<ResVK | [WeilRestrict(a * w) : w in Basis(WK), a in K]>;
RatW := RatSubspace meet ResWK;

//////////////////
// Norms 

/* assert Degree(K) eq 3; */
/* norms := {R ! Norm(iphiK(g)) : g in O1K}; */
/* W1 := sub<VK | [toVK(g) : g in norms]>; */


//////////////////
// Norms times linears

/* assert Degree(K) eq 2; */
/* norms := {R ! Norm(imp(g)) : g in O1K}; */
/* norm_times_linear := {R ! (iphi(ll) * g) : g in norms, ll in O1}; */
/* W1 := sub<V | [toV(g) : g in norm_times_linear]>; */
