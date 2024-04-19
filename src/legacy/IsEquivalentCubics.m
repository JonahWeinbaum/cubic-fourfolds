//
// NOTE: I've added a version of this function to CubicLib. Hopefully this file can now
//       be depreciated.
//

// This code determines whether 2 cubics f1, f2 are in the same PGL6-orbit, and if they are, gives
// a matrix such that f1^g = f2.

// Magma annoyingly will not compute conjugating matrix, instead we have to re-type everything
// as a permuation group using the OrbitAction intrinsic

function ConjugatingMatrix(grp, X, Y)
    map, permgp, matgp := OrbitAction(grp, {X,Y});

    bool, g := IsConjugate(permgp, map(X), map(Y)) ;
    if not IsConjugate(permgp, map(X), map(Y)) then return false, "not conjugate";
    end if;

    return bool, g @@ map;
end function;



// Now we will use our usual filtrations by irred subspaces V > Wp > W to write
// IsEquivalentCubics function

k := FiniteField(2);
G := GL(6, k);
P5<[x]> := ProjectiveSpace(k, 5);
R := CoordinateRing(P5);

V, Bit := GModule(G, R, 3);
vmap := GModuleAction(V);

W := sub<V | Bit(x[1]^3) > ;

phi := Morphism(W, V);

GV := ActionGroup(V);
GW := ActionGroup(W);

U, piU := quo<V | W>;

GU := ActionGroup(U);


//Hom from G to GL(U)
umap := GModuleAction(U);

//space of quadratic*linear, which has dimension 36
Wp := sub<V | Bit(x[1]^2*x[2])>;


Up, phip := V/Wp;
upmap := GModuleAction(Up);

GUp := ActionGroup(Up);

//Given two cubic forms f1, f2, find some matrix g in G such that (Bit(f1))^g = Bit(f2), or else return "not in same orbit". 

function IsEquivalentCubics(f1, f2)

    v1 := Bit(f1);
    vU1 := piU(v1);
    vUp1 := phip(v1);
    
    v2 := Bit(f2);
    vU2 := piU(v2);
    vUp2 := phip(v2);

    bool, gUp := ConjugatingMatrix(GUp, vUp1, vUp2);

    if bool eq false then return "not equivalent"; end if;

    vU1gUp1 := vU1^(umap(gUp @@ upmap));

    /*
      stabUp1 := Stabilizer(GUp, vUp) @@upmap;	
      Vstab1, Bitstab1 := GModule(stabUp1, R, 3);
    
    Wstab1 := sub<Vstab1 | [Vstab1! Vector(phi(w))  : w in Basis(W)]>;  
    Ustab1, pistab1 := Vstab1/Wstab1;

    ustabmap1 := GModuleAction(Ustab1);    
    GUstab1 := ActionGroup(Ustab1);
   */

    stabUp2 := Stabilizer(GUp, vUp2) @@ upmap;
    Vstab2, Bitstab2 := GModule(stabUp2, R, 3); 

    
    Wstab2 := sub<Vstab2 | [Vstab2! Vector(phi(w))  : w in Basis(W)]>;  

    Ustab2, pistab2 := Vstab2/Wstab2;

    ustabmap2 := GModuleAction(Ustab2);
    
    GUstab2 := ActionGroup(Ustab2);

    bool, gU := ConjugatingMatrix(GUstab2,  vU1gUp1, vU2);

    if bool eq false then return "not equivalent"; end if;

    stabU2 := Stabilizer(GUstab2, vU2) @@ umap;

    Vstab2, Bitstab2 := GModule(stabU2, R, 3);   

    GVstab2 := ActionGroup(Vstab2);

    bool, gV := ConjugatingMatrix(GVstab2, (v1^(vmap(gUp @@ upmap)))^(vmap(gU @@ umap)), v2);

    if bool eq false then return "not equivalent"; end if;


    return (gUp @@upmap) * (gU@@umap) * (gV@@vmap);

end function;
