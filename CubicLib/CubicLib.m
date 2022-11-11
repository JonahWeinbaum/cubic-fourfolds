//////////////////////////////////////////////////////////////////////////////////////////////////
//
// CubicLib.m
//
// Library for managing the computation of information about cubic 4-folds.
//
//////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////
//
// Constants
//
/////////////////////////////////////////////////

PATH_TO_LIB := PathToLib();

ORBIT_DATA_DIRECTORY := PATH_TO_LIB * 
                        "../../database/group_action/orbit_representatives/" *
                        "orbit_representative_in_V/";

DATA_DIRECTORY := PATH_TO_LIB * "../../database/zeta/";
CUBIC_ID_FILE := "orbrep.csv";
ISSMOOTH_FILE := "smooth.csv";
POINT_COUNTS_FILE := "point_counts.csv";
ORBIT_SIZE_FILE := "stabilizer_counts.csv";
ERROR_FILE := "error_report";

/////////////////////////////////////////////////
//
// Cubic/Bitstring conversions
//
/////////////////////////////////////////////////

CONST_k := FiniteField(2);
CONST_G := GL(6, CONST_k);
CONST_R<[x]> := PolynomialRing(CONST_k, 6);
CONST_V, CONST_Bit := GModule(CONST_G, CONST_R, 3);
    
intrinsic CubicToBitList(cubic) -> SeqEnum
{Converts a cubic into a (normalized) list of 56 bits}
    
    return Reverse([Integers() ! CONST_Bit(cubic)[i] : i in [1..56]]);
end intrinsic;

intrinsic BitListToCubic(blist) -> RngMPolElt
{Converts a list of 56 bits to a cubic}
    return blist @@ CONST_Bit;
end intrinsic;

intrinsic CubicToInt(cubic) -> RngIntElt
{}
    return Seqint(CubicToBitList(cubic), 2);
end intrinsic;

intrinsic GetBitMap() -> .
{}
    return CONST_Bit;
end intrinsic;

intrinsic StandardCubicModule() -> ModGrp
{Return CubicLib's internal GModules.}
    return CONST_V, CONST_Bit;
end intrinsic;


/////////////////////////////////////////////////
//
// Warring-like subspaces
//
/////////////////////////////////////////////////

// The Polemic subspace is the space of quadratic*linear, which has dimension 36.
// The sub<...> constructor is in the category of G-modules.
CONST_Warring := sub<CONST_V | CONST_Bit(CONST_R.1^3)> ;
CONST_Polemic := sub<CONST_V | CONST_Bit(CONST_R.1^2 * CONST_R.2)>;

CONST_Polemic_quo, CONST_to_polemic_quo := quo<CONST_V | CONST_Polemic>;
CONST_Warring_quo, CONST_to_warring_quo := quo<CONST_V | CONST_Warring>;

intrinsic WarringFiltration(v) -> ModGrpElt, ModGrpElt, ModGrpElt
{Given a vector `v` of length 56 representing a cubic, return the images
`v`, `vU`, `vUp` in the Warring filtration.

This function is very picky about its domain. See StandardCubicModule.}
    return v, CONST_to_warring_quo(v), CONST_to_polemic_quo(v);
end intrinsic;

intrinsic WarringActionMaps() -> Map, Map, Map
{}
    cub_action := GModuleAction(CONST_V);
    pol_action := GModuleAction(CONST_Polemic_quo);
    war_action := GModuleAction(CONST_Warring_quo);    
    return cub_action, pol_action, war_action;
end intrinsic;

/////////////////////////////////////////////////
//
// Orbits and stabilizers
//
/////////////////////////////////////////////////

intrinsic Stabilizer(v :: ModGrpElt) -> GrpMat
{Construct the stabilizer of `v` within the action group.}
    G := ActionGroup(Parent(v));
    return Stabilizer(G, v);    
end intrinsic;
    
intrinsic CubicStabilizer(f :: RngMPolElt) -> GrpMat
{Computes the stabilizer of a cubic form over F2 under the action of PGL_6.}

    // To efficiently compute the stabilizer, we use the standard filtration.
    // as described in the paper.
    
    // Images of f under the standard filtration.
    v, vU, vUp := WarringFiltration(CONST_Bit(f));

    // Action maps
    cub_action, pol_action, war_action := WarringActionMaps();
    
    // Compute stabilizer from the Polemic quotient.
    StabvUp := Stabilizer(vUp);
    WQres   := Restriction(CONST_Warring_quo, StabvUp @@ pol_action);

    // Compute stabilizer from the Warring quotient.
    StabvU := Stabilizer(WQres ! vU);    
    Vres   := Restriction(CONST_V, StabvU @@ war_action);

    // Compute the stabilizer.
    Stabv := Stabilizer(Vres ! v);

    // Return the result within PGL_6
    return Stabv @@ cub_action;
end intrinsic;


intrinsic CubicOrbitSize(f :: RngMPolElt) -> RngIntElt
{Computes the size of the orbit of a cubic form over F2 under the action of PGL_6.}
    Stabv := CubicStabilizer(f);
    return #CONST_G/#Stabv;    
end intrinsic;
    
/////////////////////////////////////////////////
//
// Geometry
//
/////////////////////////////////////////////////

intrinsic IsSmooth(f :: RngMPolElt) -> BoolElt
{}
    R  := Parent(f);
    PP := Proj(R);
    X := Scheme(PP, f);
    return IsNonSingular(X);
end intrinsic;


function ConjugatingMatrix(grp, X, Y)
    // This code determines whether 2 cubics f1, f2 are in the same PGL6-orbit, and if they are,
    // gives a matrix such that f1^g = f2.

    // Magma annoyingly will not compute conjugating matrix, instead we have to recast
    // everything as a permuation group using the OrbitAction intrinsic
    
    map, permgp, matgp := OrbitAction(grp, {X,Y});

    bool, g := IsConjugate(permgp, map(X), map(Y));
    
    if not IsConjugate(permgp, map(X), map(Y)) then
        return false, "not conjugate";
    end if;

    return bool, g @@ map;
end function;

intrinsic IsConjugate(G::GrpMat, v::ModGrpElt, w::ModGrpElt) -> BoolElt, GrpMatElt
{Checks if G sends v to w, and if so, returns a matrix which does so.}

    assert G cmpeq ActionGroup(Parent(v));
    
    mp, permgp, matgp := OrbitAction(G, {v,w});
    bool, g := IsConjugate(permgp, mp(v), mp(w));

    if bool then
        return bool, g @@ mp;
    else
        return false;
    end if;
end intrinsic;


intrinsic IsEquivalentCubics(f1, f2) -> BoolElt, GrpMatElt
{Given two cubic forms f1, f2, determine if there is some matrix g in G such that 
Bit(f1)^g = Bit(f2). If it exists, also return the transformation `g`.}

    // Images of f under the standard filtration.
    v1, vU1, vUp1 := WarringFiltration(CONST_Bit(f1));
    v2, vU2, vUp2 := WarringFiltration(CONST_Bit(f2));

    Wp := Parent(vUp1);
    W  := Parent(vUp2); 

    // Action maps
    cub_action, pol_action, war_action := WarringActionMaps();

    // Rename action maps so Jack's code works.
    upmap := pol_action;
    umap  := war_action;
    vmap  := cub_action;
    
    GUp := ActionGroup(Parent(vUp1));
    GU  := ActionGroup(Parent(vU1));
    
    bool, gUp := IsConjugate(GUp, vUp1, vUp2);
    if bool eq false then return false; end if;
    vU1gUp1 := vU1^(umap(gUp @@ pol_action));

    // This is v1 moved over the point [v2] in the Polemic quotient
        
    // Compute stabilizer from the Polemic quotient.
    StabvUp := Stabilizer(vUp2);
    WQres   := Restriction(CONST_Warring_quo, StabvUp @@ pol_action);
    
    bool, gU := IsConjugate(ActionGroup(WQres), vU1gUp1, vU2);
    if bool eq false then return false; end if;
    next_thing := (v1^(vmap(gUp @@ upmap)))^(vmap(gU @@ umap));

    
    // Compute the final refinement, if it exists.
    StabvU2 := Stabilizer(WQres ! vU2);
    Vres    := Restriction(CONST_V, StabvU2 @@ war_action);
    
    bool, gV := IsConjugate(ActionGroup(Vres), next_thing, v2);
    if bool eq false then return false; end if;

    return true, (gUp @@ upmap) * (gU @@ umap) * (gV @@ vmap);
end intrinsic;



/////////////////////////////////////////////////
//
// Group theory
//
/////////////////////////////////////////////////

// Count Orbits using Burnside.
intrinsic CountOrbits(Gmod) -> RngIntElt
{Counts the orbits of a group acting on a G-module via Burnside's lemma.}
    rho := Representation(Gmod);
    G := Group(Gmod);
    return &+[#Eigenspace(rho(c[3]),1) * c[2] : c in ConjugacyClasses(G) ]/#G;
end intrinsic;


/////////////////////////////////////////////////
//
// Data Processing
//
/////////////////////////////////////////////////

intrinsic serialize(f :: RngMPolElt) -> .
{Takes a cubic and returns a customized serialization into bytes.}
    byteseq := [];
    for i in [0..6] do
        n := 1;
        bitstring := 0;
        for j in [1..8] do 
            if f[i*8 + j] eq 1 then
                bitstring := BitwiseOr(bitstring, ShiftLeft(1, 8-n));
            end if;
            n := n + 1;
        end for;
        Append(~byteseq, bitstring);
    end for;

    return byteseq;
end intrinsic;

intrinsic deserialize(byteseq) -> RngMPolElt
{Takes a byte sequence and returns the corresponding cubic polynomial.}
    seq := [];
    for byte in byteseq do
        seq cat:= [ BitwiseAnd(ShiftRight(byte, 8-i), 1): i in [1..8]];
    end for;
    
    return seq;
end intrinsic;

intrinsic LoadCubicOrbitData(: RemoveZero:=true, Flat:=false, Quick:=false, BitList:=false)
          -> SeqEnum
{Loads the precomputed orbit data.}

    ZEROCUBIC := [0 : i in [1..56]];
    print "loading data...";

    if Quick then
	range := [1..2];
    else
	range := [1..85];
    end if;
    
    orbdata := [];
    time for k in range do
	     name2 := Sprintf(ORBIT_DATA_DIRECTORY * "orbitreps-%o.data", k);
	     seq := [];
	     o2 := []; 
	     file := Read(name2);

	     currloc := 1;
	     while true do
		 //Deserialize 7 Bytes and Advance File Pointer   
		 for i in [1..7] do
		     byte := StringToCode(file[currloc]);
		     
		     seq cat:= [byte];   
		     
		     //Break at EOF
		     if (currloc eq #file) then
			 break;
		     end if;
		     currloc := currloc + 1;
		 end for;

		 //Deserialize
		 x := deserialize(seq);
		 if not (RemoveZero and x eq ZEROCUBIC) then
		     Append(~o2, x);
		 end if;		 
		 seq := [];

		 //Break at EOF
		 if (currloc eq #file) then
		     break;
		 end if;
	     end while;

	     orbdata := orbdata cat [o2];
	     k;
	 end for;

    if not BitList then
        orbdata := [[BitListToCubic(f) : f in lst] : lst in orbdata];
    end if;
    
    if Flat then
	return &cat orbdata;
    else
	return orbdata;
    end if;
end intrinsic;


/////////////////////////////////////////////////
//
// Database
//
/////////////////////////////////////////////////

// Write stuff to file.

// We would like to include records of the following form:
//
// 1. The orbit representative a, b such that it is orbdata[a,b];
// 2. Whether the orbit representative is smooth.
// 3. The sequence of point counts.

// 4. Set of lines on it.
// 5. Set of planes on it.

intrinsic WriteZetaData(i, issmooth, pointcounts) -> RngIntElt
{}
    Write(DATA_DIRECTORY * CUBIC_ID_FILE, Sprintf("%o", i));
    Write(DATA_DIRECTORY * ISSMOOTH_FILE, Sprintf("%o, %o", i, issmooth));
    Write(DATA_DIRECTORY * POINT_COUNTS_FILE, Sprintf("%o, %o", i, pointcounts));    
    return 0;
end intrinsic;

intrinsic WriteOrbitSizeData(i, stabilizersize) -> RngIntElt
{}
    Write(DATA_DIRECTORY * ORBIT_SIZE_FILE, Sprintf("%o, %o", i, stabilizersize));
    return 0;
end intrinsic;

intrinsic ReportError(index, err)
{Write a report of the error to the file.}
    Write(DATA_DIRECTORY * ERROR_FILE, Sprintf("%o, %o", index, err));
    return;
end intrinsic;

/////////////////////////////////////////////////
//
// C++ interaction.
//
/////////////////////////////////////////////////

intrinsic CppInputString(h) -> MonStgElt
{Given a quantity `h`, return the string needed for Nick Addington's C++ program.}
    if h eq 0 then return "0";
    end if;
    ret := "";
    for m in Monomials(h) do
        j := [Degree(m, Parent(h).i) : i in [1..Rank(Parent(h))] ];
        str := "";
        for k in [0 .. #j-1] do
            for l in [1..j[k+1]] do
                if str eq "" then str := "y_" cat Sprint(k);
                else str := "mult[" cat str cat "][" cat "y_" cat Sprint(k) cat "]";
                end if;
            end for;
            
        end for;
        
        if ret eq "" then ret := str;
        else ret := ret cat " ^ " cat str;
        end if;
    end for;
    
    return ret;
end intrinsic;


/////////////////////////////////////////////////
//
// find the characteristic polynomial of frobenius on nonprimitive cohomology.
// input list of point counts over F_2^k, k = 1,..., 11.
/////////////////////////////////////////////////

CONSTQt_<t> := PolynomialRing(Rationals());
CONSTCs_<s> := PolynomialRing(ComplexField(30));

intrinsic Charpoly(list::SeqEnum[RngIntElt]) -> RngUPolElt, RngIntElt, MonStgElt
{input list of point counts over F_2^k, k = 1,..., 11, 
returns the characteristic polynomial of Frobenius acting on nonprimitive cohomology.}

    tr := [(list[m] - 1 - 2^m - 4^m - 8^m - 16^m)/4^m : m in [1..11]];
    c := [];
    c[1] := -tr[1];
    for k in [2..11] do
        c[k] := (tr[k] + &+[c[i]*tr[k-i] : i in [1..k-1]] )/(-k);
    end for;
    
    p := t^22 + &+[c[i]*t^(22-i) : i in [1..11] ];
    g := CONSTQt_ ! (t^22 * Evaluate(p, 1/t));
    
    if c[11] ne 0 then
        ret := p + Modexp(g, 1, t^11);
        return ret, 1;
    else
        C := ComplexField(30);
        poly1 := p + Modexp(g, 1, t^11);
        roots1 := Roots(Evaluate(poly1, s));
        mods1 := [Modulus(roots1[i][1]) : i in [1..#roots1]];
        
        if &and [Abs(mods1[i] - 1) lt 1E-30 : i in [1..#mods1] ] then
            return poly1, 1, "middle coeff zero";
        else
            poly2 := p - Modexp(g, 1, t^11);
            roots2 := Roots(Evaluate(poly2, s));
            mods2 := [Modulus(roots2[i][1]) : i in [1..#roots2]];
            
            if &and [Abs(mods2[i] - 1) lt 1E-30 : i in [1..#mods2] ] then
                return poly2, -1, "middle coeff zero";
            else
                error "Something went wrong in zeta function computation.";
            end if;
        end if;

    end if;
    
end intrinsic;
    
    
