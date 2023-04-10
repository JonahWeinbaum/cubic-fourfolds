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

DATABASE_DIRECTORY := PATH_TO_LIB * "../../database/";
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

intrinsic StandardCubicRing() -> ModGrp
{Return CubicLib's internal GModules.}
    return CONST_R;
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


intrinsic IsConjugate(G::GrpMat, v::ModGrpElt, w::ModGrpElt) -> BoolElt, GrpMatElt
{Checks if G sends v to w, and if so, returns a matrix which does so.}

    mp, permgp, matgp := OrbitAction(G, {v,w});
    bool, g := IsConjugate(permgp, mp(v), mp(w));

    if bool then
        return bool, g @@ mp;
    else
      return false, _;
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

CONST_Max_Feasible_Orbits := 10^14;

// Count Orbits using Burnside.
intrinsic CountOrbits(Gmod) -> RngIntElt
{Counts the orbits of a group acting on a G-module via Burnside's lemma.}
    rho := Representation(Gmod);
    G := Group(Gmod);
    return &+[#Eigenspace(rho(c[3]),1) * c[2] : c in ConjugacyClasses(G) ]/#G;
end intrinsic;

intrinsic GLModule(n::RngIntElt, d::RngIntElt, q::RngIntElt) -> ModGrp
{Returns the module for the action of GL(n+1, q) on spaces of polynomials of degree d.}
    G := GL(n+1, q);
    R := PolynomialRing(GF(q), n+1);
    V := GModule(G, R, d);
    return V;
end intrinsic;

intrinsic IsFeasible(n::RngIntElt, d::RngIntElt, q::RngIntElt
                     : Ncores := 100,
                       TimeLimit := 525600 * 60) -> BoolElt
{A Heuristic to determines whether orbit enumeration is feasible for forms
of degree d in P^n over Fq. Our heuristic is based on current techniques 
based on reasonable hardware from 2022.

NOTE: Sometimes this function will produce a false negative, because it 
doesn't find the optimal composition series.
}
    G := GL(n+1, q);
    R := PolynomialRing(GF(q), n+1);
    V := GModule(G, R, d);    
    return IsFeasible(V : Ncores:=Ncores, TimeLimit:=TimeLimit);
end intrinsic;

intrinsic IsFeasible(V::ModGrp
                     : Ncores := 100,
                       TimeLimit := 525600 * 60) -> BoolElt
{}
    if CountOrbits(V) gt CONST_Max_Feasible_Orbits then
        return false;
    end if;
    
    // This can take a while for really big modules. The answer is probably no
    // in this case anyways.
    found_series := CandidateCompositionSeries(V);
    return &or [IsFeasibleUnionFind(V) : V in {comp[1] : comp in found_series}];
end intrinsic;

intrinsic CandidateCompositionSeries(V::ModGrp) -> SeqEnum
{}
    // Look for 10 composition series.
    found_series := [];
    found_dims := [];
    
    wait_count := 0;    
    while wait_count lt 10 do
        comp := CompositionSeries(V);
        dims := [Dimension(W) : W in comp];
        
        if not dims in found_dims then
            Append(~found_dims, dims);
            Append(~found_series, comp);
            wait_count := 0;
        else
            wait_count +:= 1;
        end if;
    end while;

    return found_series;
end intrinsic;

intrinsic IsFeasibleUnionFind(n::RngIntElt, d::RngIntElt, q::RngIntElt
                              : Ncores := 100,
                                TimeLimit := 525600 * 60) -> BoolElt
{A Heuristic to determines whether orbit enumeration is feasible for forms
of degree d in P^n over Fq. Our heuristic is based on using union-find 
based on reasonable hardware from 2022.}

    if n gt 2 and d gt 10 then
        return false;
    end if;

    G := GL(n+1, q);
    R := PolynomialRing(GF(q), n+1);    
    V := GModule(G, R, d);

    return IsFeasibleUnionFind(V : Ncores:=Ncores, TimeLimit:=TimeLimit);
end intrinsic;

intrinsic IsFeasibleUnionFind(V::ModGrp
                              : Ncores := 100,
                                TimeLimit := 525600 * 60) -> BoolElt
{}
    // Note the number of seconds in a year is 525600 * 60.

    if CountOrbits(V) gt CONST_Max_Feasible_Orbits then
        return false;
    end if;
    
    // If the dimension is larger than 60, the answer is no.
    if Dimension(V) gt 60 then return false; end if;
    
    // We use the complexity estimate that union find is O(N).
    // The group PGL_{n+1} can be assumed to act (generically) freely, since
    // cubics with larger stabilizers are generally sparse within the family.
    
    // Time matrix-vector multiplication.
    G := Group(V);
    gens := Generators(G);
    ntrials := 10^4;
    
    t1 := Cputime();
    for i in [1..ntrials] do
        v := Random(V);
        g := Random(gens);
        gv := v * g;
    end for;
    t2 := Cputime();
    avg_step_time := (t2-t1)/ntrials;

    // We now compare to hardware availability. We assume that time complexity
    // is the main contraint.
    q := #BaseRing(V);
    N := ExactQuotient(#V-1, q-1);

    return N * avg_step_time / Ncores lt TimeLimit;
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

intrinsic LoadCubicOrbitData(: RemoveZero:=true, Flat:=false, Quick:=false, BitList:=false,
                               OnlySmooth:=false, Verbose:=false)
          -> SeqEnum
{Loads the precomputed orbit data.}

    ZEROCUBIC := [0 : i in [1..56]];
    if Verbose then print "loading data..."; end if;

    if Quick then
	range := [1..2];
    else
	range := [1..85];
    end if;
    
    orbdata := [];
    for k in range do
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
	if Verbose then print k; end if;
    end for;

    if not BitList then
        orbdata := [[BitListToCubic(f) : f in lst] : lst in orbdata];

        if OnlySmooth then
            orbdata := [[f : f in lst | IsSmooth(f)] : lst in orbdata];
        end if;
    else
        if OnlySmooth then
            orbdata := [[b : b in lst | IsSmooth(BitListToCubic(b))] : lst in orbdata];
        end if;
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

intrinsic DatabaseDirectory() -> MonStgElt
{Return the absolute path to the database folder.}
    return DATABASE_DIRECTORY;
end intrinsic;

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

intrinsic ReadCSV(fname : DataPath:=DatabaseDirectory()) -> Assoc
{}
    io := Open(DataPath * fname, "r");
    data := AssociativeArray();
    while true do
        line := Gets(io);
        if IsEof(line) then break; end if;
        args := eval("[* " * line * " *]");
        i := args[1];
        if #args eq 1 then
            data[i] := true; // Just mark the key.
        elif #args eq 2 then
            data[i] := args[2];
        else
            data[i] := args[2..#args];
        end if;
    end while;
    return data;
end intrinsic;

intrinsic Index(A::Assoc, x::Any) -> Any
{}
    return {k : k in Keys(A) | A[k] eq x};
end intrinsic;

/////////////////////////////////////////////////
//
// C++ interaction.
//
/////////////////////////////////////////////////

intrinsic CppInputString(h) -> MonStgElt
{Given a quantity `h`, return the string needed for Nick Addington's C++ program.}
    if h eq 0 then return "0"; end if;
    ret := "";
    for m in Monomials(h) do
        j := [Degree(m, Parent(h).i) : i in [1..Rank(Parent(h))] ];
        str := "";
        for k in [0 .. #j-1] do
            for l in [1..j[k+1]] do
                if str eq "" then
                    str := Sprintf("y_%o", k);
                else
                    str := Sprintf("ff2k_mult(%o, y_%o)", str, k);
                    // str := "mult[" cat str cat "][" cat "y_" cat Sprint(k) cat "]";
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
// Macaulay2 interaction
//
/////////////////////////////////////////////////

intrinsic FunctionalEquationSign(f) -> RngIntElt
{Given a cubic fourfold `f` over F2, compute the sign of the functional equation.
This function requires that either Macaulay 2 or Sage is installed.}

    // Testing functionality with the Sage fallback first.
    R<x,y,z,u,v,w> := PolynomialRing(BaseRing(Parent(f)), 6);
    fR := R ! f;

    fstr := Sprint(fR);
    try
	// Write to input file and execute M2.
	Write("macaulay2/input.m2", fstr : Overwrite);
	cmd := Sprint("M2 --script macaulay2/compute_sign.m2");
        disc := eval Read(POpen(cmd, "r"));
        
    catch e
        print "Warning: Using Sage fallback. Please ensure Macaulay2 can be launched with `M2`.";
        // If we use sage, there is a missing factor of
        //     3^21 = 3 mod 8.
        // compared to the Macaulay implementation.        
        cmd := Sprintf("sage -c 'f=\"%o\"; load(\"macaulay2/compute_sign.sage\")'", fstr);
        disc := 3 * eval Read(POpen(cmd, "r"));
    end try;
    
    // Once the discriminant is found, compute the sign.
    if disc mod 8 eq 3 then
        return -1;
    elif disc mod 8 eq 1 then
        return 1;
    else
        error "Unexpected discriminant residue.";
    end if;
    
end intrinsic;

/////////////////////////////////////////////////
//
// find the characteristic polynomial of frobenius on nonprimitive cohomology.
// input list of point counts over F_2^k, k = 1,..., 11.
//
/////////////////////////////////////////////////

CONSTQt_<t> := PolynomialRing(Rationals());
CONSTCs_<s> := PolynomialRing(ComplexField(30));

function WeilPolynomialFromHalf(coeffs, sign)
    head := coeffs[0+1..11+1];
    tail := [sign*x : x in Reverse(coeffs[0+1..10+1])];
    return Polynomial(Reverse(head cat tail)); // I'm not a fan...constant term should be 1.
end function;

intrinsic BothWeilPolynomials(coeffs) -> Tup
{}
    a := WeilPolynomialFromHalf(coeffs,  1);
    b := WeilPolynomialFromHalf(coeffs, -1);
    return <a,b>;
end intrinsic;


intrinsic HalfWeil(counts::SeqEnum[RngIntElt]) -> SeqEnum[FldRatElt]
{}
    // Divide by the zeta function for P4. (Recall there is a Log.)
    // Also, dividing by 4^m amounts to taking a Tate twist (evaluating t => t/4).
    tr := [(counts[m] - 1 - 2^m - 4^m - 8^m - 16^m)/4^m : m in [1..11]];

    // Apply the exponential map. Note that c[i] is the i-th coefficient.
    c := [Rationals() ! -2^31 : i in [1..11]]; // Initialize
    c[1] := -tr[1];
    for k in [2..11] do
        c[k] := (tr[k] + &+[c[i]*tr[k-i] : i in [1..k-1]])/(-k);
    end for;

    // Append the initial 1. (Note: this might be a sign different than the old code.)
    return [1] cat c;
end intrinsic;

intrinsic IsSignAmbiguous(halfWeil) -> BoolElt, RngIntElt
{Determine if there is an ambiguity in the functional equation. Also return the
"Badness" rating -- the largest integer such that the 11+ith coefficient is zero.}

    if halfWeil[12] ne 0 then return false, -1; end if;
    
    // Try both signs. See if only one works.
    poly1 := WeilPolynomialFromHalf(halfWeil,  1);
    poly2 := WeilPolynomialFromHalf(halfWeil, -1);

    isWeil1 := HasAllRootsOnUnitCircle(poly1);
    isWeil2 := HasAllRootsOnUnitCircle(poly2);

    if isWeil1 and isWeil2 then
        badness := 0;
        while true do
            if halfWeil[12-badness] ne 0 then break; end if;
            badness +:=1;
        end while;
        return true, badness-1;
    else
        return false, -1;
    end if;
end intrinsic;

intrinsic Badness(halfWeil) -> RngIntElt
{}
    a, b := IsSignAmbiguous(halfWeil);
    return b;
end intrinsic;

intrinsic DoesArtinTateHelp(halfWeil) -> BoolElt
{Check if the value at -1 is a square or twice a square.}
    atv := Integers() ! (2*Evaluate(Polynomial(halfWeil), -1));
    return not (IsSquare(atv) or IsSquare(2*atv));
end intrinsic;

intrinsic Charpoly(list::SeqEnum[RngIntElt] : FailIfAmbiguous:=true) -> RngUPolElt
{input list of point counts over F_2^k, k = 1,..., 11, 
returns the characteristic polynomial of Frobenius acting on nonprimitive cohomology.}
    
    // p := t^22 + &+[c[i]*t^(22-i) : i in [1..11] ];
    // g := CONSTQt_ ! (t^22 * Evaluate(p, 1/t));

    c := HalfWeil(list);
    
    if c[11+1] ne 0 then
        return WeilPolynomialFromHalf(c, 1);
    end if;


    // Try both signs. See if only one works.
    poly1 := WeilPolynomialFromHalf(c,  1);
    poly2 := WeilPolynomialFromHalf(c, -1);

    isWeil1 := HasAllRootsOnUnitCircle(poly1);
    isWeil2 := HasAllRootsOnUnitCircle(poly2);
    
    if isWeil1 and not isWeil2 then
        return poly1;
    elif isWeil2 and not isWeil1 then
        return poly2;

    elif isWeil1 and isWeil2 then

        // TODO: We should have a file keeping track of the higher point counts,
        // for the cubics that need them.
        
        // Note: should include an optimization regarding which point
        // count to request.

        if FailIfAmbiguous then
            error "Ambiguous sign of functional equation. Cannot determine zeta functions.";
        else
            return -1; // Obviously wrong.
        end if;


        error "Not implemented.";
        
        // TODO: The indexing is wrong
        /*
          p12  := PointCounts(cubic, 12);
          tr12 := (p12 - 1 - 2^12 - 4^12 - 8^12 - 16^12)/4^12;
          c12  := (tr12 + &+[c[i]*tr[12-i] : i in [1..12-1]])/(-12);

        if Coefficient(poly1, 12) eq c12 and c12 ne 0 then 
        return poly1; 
        elif Coefficient(poly2, 12) eq c12 and c12 ne 0 then 
        return poly2;
        else
        error "counting up to 2^12 does not cut it";
        end if;
       */
    else
        error "No valid sign choice for given point counts.";
    end if;    
end intrinsic;


/////////////////////////////////////////////////
//
// Weil polynomial manipulations.
//
/////////////////////////////////////////////////

// 100 is a sufficient cut-off point to list all cyclotomic polynomials of degree at most 24.
CONST_cyclo := [f : f in [CyclotomicPolynomial(i) : i in [1..100]] | Degree(f) le 24];

intrinsic IsWeilPolynomial(f::RngUPolElt : Proof:=false) -> BoolElt
{Checks if `f` is a Weil polynomial.}
    return HasAllRootsOnUnitCircle(f);
end intrinsic;
    

intrinsic TranscendentalFactor(f :: RngUPolElt) -> RngUPolElt
{Given a Weil polynomial, return the transcendental factor. That is,
remove all cyclotomic factors from f.}
    
    if f eq 0 then return 0; end if;
    q := f;
    for g in CONST_cyclo do
        while true do
            nextq, r := Quotrem(q, g);
            if r ne 0 then break; end if;
            q := nextq;
        end while;
    end for;
    return q;
end intrinsic;

intrinsic TranscendentalRank(f :: RngUPolElt) -> RngUPolElt
{Given a Weil polynomial, return the rank of the transcendental part.}
    return Degree(TranscendentalFactor(f));
end intrinsic;


intrinsic IrrationalFactor(f :: RngUPolElt) -> RngUPolElt
{Given a Weil polynomial remove all factors of 1-T from f.}
    R<t> := Parent(f);
    if f eq 0 then return 0; end if;
    q := f;
    while true do
        nextq, r := Quotrem(q, 1-t);
        if r ne 0 then break; end if;
        q := nextq;
    end while;
    return q;    
end intrinsic;

intrinsic TateTwist(f::RngUPolElt, j::FldRatElt : q:=2) -> RngUPolElt
{Scale the roots of f by 2^(-j)}
    require Denominator(j) eq 1 : "Not implmented for non-integral j.";
    
    _<t> := Parent(f);
    return Evaluate(f, q^(-j) * t);
end intrinsic;

intrinsic TateTwist(f::RngUPolElt, j::RngIntElt : q:=2) -> RngUPolElt
{Scale the roots of f by 2^(-j)}
    _<t> := Parent(f);
    return Evaluate(f, q^(-j) * t);
end intrinsic;

intrinsic CubicWeilPolynomialToLPolynomial(f::RngUPolElt : q:=2) -> RngUPolElt
{Converts the Weil polynomial to the zeta function.}
g := TateTwist(f, -2 : q:=2); // Untwist.
    return Coefficient(g,0) eq -1 select -g else g;
end intrinsic;

intrinsic CWTL(f) -> Any
{Alias for CubicWeilPolynomialToLPolynomial.}
    return CubicWeilPolynomialToLPolynomial(f);
end intrinsic;

intrinsic CubicLPolynomialToPointCounts(f::RngUPolElt) -> SeqEnum
{Given the L-polynomial, return the list of 22 point counts over Fq, for
q = 1, ..., 22.}
    PP<t> := PowerSeriesRing(Rationals(), 30);
    LP4inv := (1-t) * (1-2*t) * (1-4*t) * (1-8*t) * (1-16*t);
    F := PP ! f;
    g := -Log(F * LP4inv);
    coeffs := Eltseq(g)[1..22];
    return [Integers() ! (n*coeffs[n]) : n in [1..22]];
end intrinsic;

intrinsic CubicWeilPolynomialToPointCounts(f::RngUPolElt : q:=2) -> SeqEnum
{}
    return CubicLPolynomialToPointCounts(CubicWeilPolynomialToLPolynomial(f : q:=q));
end intrinsic;

intrinsic ArtinTateValue(f::RngUPolElt : Weil:=false) -> FldRatElt
{Given the L-polynomial of Frobenius acting on H^4(X, QQ_l), return the
Artin-tate value. If the parameter `Weil:=true`, then f is assumed to
be a Weil polynomial and is normalized.}

    g := Weil select CubicWeilPolynomialToLPolynomial(f) else f;

    R<t> := Parent(g);
    
    // See the paper of Komodo.
    lvalue := Evaluate(g/(1-4*t)^Valuation(g, 1-4*t), 1/4);

    // This is the special value of the zeta function ** divided by q^chi(X, OX, 2) **
    return 1/16 * 8/9 * 1/lvalue;

    // The numerator is always 1, which is a good sign for sure.

    // Our guess is off by a factor of 9. Not sure what that's about.

    // H^{3,1} might collapse (by 1) in the supersingular case.

end intrinsic;


intrinsic LoadKSDatabase(: Normalize:=true) -> SetEnum
{}
    fname := DATABASE_DIRECTORY * "zeta_functions/k3f2-lines.txt";
    
    F := Open(fname, "r");
    kspolys := {};
    
    while true do
        l := Gets(F);
        if IsEof(l) then break; end if;

        lst := eval(l);
        if Normalize then
            Include(~kspolys, Polynomial(Rationals(), lst)/2);
        else
            Include(~kspolys, Polynomial(lst));
        end if;
    end while;

    return kspolys;
end intrinsic;


intrinsic IsOrdinary(f :: RngUPolElt : Weil:=false) -> Any
{WARNING: This function is not implemented correctly.}    
    g := Weil select CubicWeilPolynomialToLPolynomial(f) else f;
    NP := LowerVertices(NewtonPolygon(g, 2));

    hodgepolygon := [ <0, 0>, <1, 1>, <21, 41>, <22, 44>];
    return NP eq hodgepolygon;
end intrinsic;
