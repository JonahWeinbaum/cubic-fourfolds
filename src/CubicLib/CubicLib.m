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

DATABASE_DIRECTORY := PATH_TO_LIB * "../../database/";
ZETA_DIRECTORY := DATABASE_DIRECTORY * "zeta_functions/";
GROUP_ACTION_DIRECTORY := DATABASE_DIRECTORY * "group_action/";
ORBIT_DATA_DIRECTORY := GROUP_ACTION_DIRECTORY * 
			"orbit_representatives/" *
                        "orbit_representative_in_V/";

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
                     : CheckOrbits:=true,
                       Ncores := 100,
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
    return IsFeasible(V : CheckOrbits:=CheckOrbits, Ncores:=Ncores, TimeLimit:=TimeLimit);
end intrinsic;

function IsFeasibleFactors(factors)
    return &and [IsFeasibleUnionFind(fac : CheckOrbits:=false) : fac in factors];
end function;

intrinsic IsFeasible(V::ModGrp                            
                     : CheckOrbits := true,
                       Ncores := 100,
                       TimeLimit := 525600 * 60) -> BoolElt
{}
    if CheckOrbits and CountOrbits(V) gt CONST_Max_Feasible_Orbits then
        return false, "Too many orbits";
    end if;
    
    // This can take a while for really big modules. The answer is probably no
    // in this case anyways.
    found_factors := CandidateCompositionFactors(V);
    return &or [IsFeasibleFactors(factors) : factors in found_factors];
end intrinsic;

intrinsic CandidateCompositionFactors(V::ModGrp) -> SeqEnum
{}
    // Look for 10 composition series.
    found_factors := [];
    found_dims := [];
    
    wait_count := 0;    
    while wait_count lt 10 do
        comp, factors := CompositionSeries(V);
        dims := [Dimension(W) : W in comp];
        
        if not dims in found_dims then
            Append(~found_dims, dims);
            Append(~found_factors, factors);
            wait_count := 0;
        else
            wait_count +:= 1;
        end if;
    end while;

    return found_factors;
end intrinsic;

intrinsic IsFeasibleUnionFind(n::RngIntElt, d::RngIntElt, q::RngIntElt
                              : CheckOrbits:=true,
                                Ncores := 100,
                                TimeLimit := 525600 * 60) -> BoolElt
{A Heuristic to determines whether orbit enumeration is feasible for forms
of degree d in P^n over Fq. Our heuristic is based on using union-find 
based on reasonable hardware from 2022.}

    if n gt 2 and d gt 10 then
        return false, "Too many orbits";
    end if;

    G := GL(n+1, q);
    R := PolynomialRing(GF(q), n+1);    
    V := GModule(G, R, d);

    return IsFeasibleUnionFind(V : CheckOrbits:=CheckOrbits,
                                   Ncores:=Ncores, TimeLimit:=TimeLimit);
end intrinsic;

intrinsic IsFeasibleUnionFind(V::ModGrp
                              : CheckOrbits:=true,
                                Ncores := 100,
                                TimeLimit := 525600 * 60) -> BoolElt
{}
    // Note the number of seconds in a year is 525600 * 60.

    if CheckOrbits and CountOrbits(V) gt CONST_Max_Feasible_Orbits then
        return false, "Too many orbits";
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

intrinsic serialize(f :: ModGrpElt) -> .
{Given a module element representing a cubic over F2, returns a customized serialization into bytes.}
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

intrinsic deserialize(byteseq) -> .
{Takes a byte sequence and returns the corresponding cubic polynomial.}
    seq := [];
    for byte in byteseq do
        seq cat:= [ BitwiseAnd(ShiftRight(byte, 8-i), 1): i in [1..8]];
    end for;
    
    return seq;
end intrinsic;


/////////////////////////////////////////////////
//
// Database
//
///////////////////////////////////////////////////

intrinsic _ManageBitListSmoothOptions(orbdata : BitList:=false, OnlySmooth:=false) -> SeqEnum
{}
    if BitList and not OnlySmooth then
        return orbdata;
    end if;

    // In all other cases we need the cubics.
    orbdata_cubics := [[BitListToCubic(b) : b in lst] : lst in orbdata];

    // Setup smoothness filterer based on options
    if OnlySmooth then
        // Check for the cache. Otherwise, use IsSmooth.
        try
            fname := "differentiability/smooth/smooth.csv";
            smoothIndices := Keys(ReadCSVWithTypes(fname, RngIntElt, MonStgElt));
            filterer := func<i, j, seen | (seen + j) in smoothIndices>;
        
        catch e
            filterer := func<i, j, seen | IsSmooth(orbdata_cubics[i][j])>;
        end try;
        
    else
        filterer := func<i, j, seen | true>;
    end if;

    // Select return type via bitlist option
    if BitList then
        selector := orbdata;
    else
        selector := orbdata_cubics;
    end if;

    // Iterate to construct the new object.
    seenElements := 0;
    orbdataNew := [];
    for i in [1..#orbdata] do
        lstNew := [];
        for j in [1..#orbdata[i]] do
            // Get the right things associated to (i,j)
            if filterer(i,j,seenElements) then
                Append(~lstNew, selector[i][j]);
            end if;                       
        end for;

        seenElements +:= #orbdata[i];
        Append(~orbdataNew, lstNew);
    end for;

    // Return
    return orbdataNew;
end intrinsic;


intrinsic LoadCubicOrbitData(: RemoveZero:=true, Flat:=false, Quick:=false, BitList:=false,
                               OnlySmooth:=false, Verbose:=false)
          -> SeqEnum
{Loads the precomputed orbit data.

Parameters:

    RemoveZero -- Whether or not to keep the zero orbit.
    Flat       -- If true, returns a single flattened list of orbit representatives.
                  Otherwise, a list of lists is returned, where each list corresponds to an orbit
                  over one of the cosets from the filtration.

    Quick      -- If true, return only the initial bit of the data. Useful if one does not want to
                  wait for the entire dataset to load.

    BitList    -- If true, returns the cubics as a sequence of coefficients. Otherwise, returns
                  cubic polynomials.

    OnlySmooth -- If true, returns only the smooth cubics.
    Verbose    -- Verbose parameter.
}

    ZEROCUBIC := [0 : i in [1..56]];
    if Verbose then print "loading data..."; end if;
    
    if Quick then range := [1..2]; else range := [1..85]; end if;
    
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

    // Filter by smoothness and apply Bitlist.
    orbdata := _ManageBitListSmoothOptions(orbdata : BitList:=BitList, OnlySmooth:=OnlySmooth);
        
    if Flat then
	return &cat orbdata;
    else
	return orbdata;
    end if;
end intrinsic;


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

intrinsic OrbitDataDirectory() -> MonStgElt
{Return the absolute path to the orbit data folder.}
    return ORBIT_DATA_DIRECTORY;
end intrinsic;

intrinsic WriteZetaData(i, issmooth, pointcounts) -> RngIntElt
{}
    Write(ZETA_DIRECTORY * CUBIC_ID_FILE, Sprintf("%o", i));
    Write(ZETA_DIRECTORY * ISSMOOTH_FILE, Sprintf("%o, %o", i, issmooth));
    Write(ZETA_DIRECTORY * POINT_COUNTS_FILE, Sprintf("%o, %o", i, pointcounts));    
    return 0;
end intrinsic;

intrinsic ReadZetaFunctions() -> Assoc
{Load the list of Zeta functions from the database.}
    fname := "zeta_functions/zeta_coefficients.csv";
    R<t> := PolynomialRing(Rationals());

    // zetafunctions[k] will be Q(t), the characterisitc polynomial of frobenius
    // acting on on H^4_prim(X, Q_l(2)).
    A := ReadCSV(fname);
    B := AssociativeArray();
    for k in Keys(A) do
        B[k] := R ! Polynomial(A[k]);
    end for;
    return B;
end intrinsic;


intrinsic WriteOrbitSizeData(i, stabilizersize) -> RngIntElt
{}
    Write(GROUP_ACTION_DIRECTORY * ORBIT_SIZE_FILE, Sprintf("%o, %o", i, stabilizersize));
    return 0;
end intrinsic;

intrinsic ReadOrbitSizeData() -> Assoc
{}
    fname := "group_action/stabilizers_info/stabilizer_counts.csv";
    A := ReadCSV(fname);
    return A;
end intrinsic;

intrinsic ReadStabilizerSizeData() -> Assoc
{}
    n := #GL(6, FiniteField(2));
    f := func<size | ExactQuotient(n, size)>;
    fname := "group_action/stabilizers_info/stabilizer_counts.csv";
    return ReadCSV(fname : FunctionOnLoad:=f);
end intrinsic;

intrinsic ReportError(index, err)
{Write a report of the error to the file.}
    Write(ZETA_DIRECTORY * ERROR_FILE, Sprintf("%o, %o", index, err));
    return;
end intrinsic;

intrinsic DetermineReadTypes(DataPath, fname) -> Tup
{}
    io := Open(DataPath * fname, "r");
    data := AssociativeArray();
    
    line := Gets(io);

    // If the file is empty, the return value is moot.
    if IsEof(line) then
        return MonStgElt;
    end if;
    
    args := eval("[* " * line * " *]");

    return <Type(x) : x in args>;
end intrinsic;

intrinsic StringToSeqEnum(str) -> SeqEnum
{}
    // Strip off the left/right brackets.
    leftBrac  := Index(str, "[");
    rightBrac := Index(str, "]");

    // Parse the integer sequence.
    lst := Split(str[leftBrac+1..rightBrac-1], ",");    
    return [StringToRational(x) : x in lst];
end intrinsic;

intrinsic SplitAtFirstComma(str) -> MonStgElt, MonStgElt
{}
    i := Index(str, ",");
    keystr   := str[1..i-1];
    valuestr := str[i+1..#str];
    return keystr, valuestr;
end intrinsic;

intrinsic GenericSplitterParser(str) -> Any, Any
{}
    args := eval("[* " * str * " *]");
    assert #args eq 2;
    return args[1], args[2];
end intrinsic;

intrinsic GenericParser(str) -> Any
{}
    return eval str;
end intrinsic;
    
intrinsic ReadCSVWithTypes(fname, keyType::Cat, valueType::Cat :
                           DataPath:=DatabaseDirectory(),
                           FunctionOnLoad := func<x|x>) -> Assoc
{}
    io := Open(DataPath * fname, "r");
    data := AssociativeArray();

    case keyType:
    when RngIntElt: KeyParser := StringToInteger; Splitter:=SplitAtFirstComma;
    else: KeyParser := func<x|x>; Splitter:=GenericSplitterParser;
    end case;
                
    case valueType:
    when RngIntElt: ValueParser := StringToRational;
    when FldRatElt: ValueParser := StringToRational;
    when SeqEnum:   ValueParser := StringToSeqEnum;
    when MonStgElt: ValueParser := func<x|x>;
    else: ValueParser := GenericParser;
    end case;

    // Override the Value parser if we have a generic Keys
    if Splitter eq GenericSplitterParser then
        ValueParser := func<x|x>;
    end if;

    while true do
        line := Gets(io);
        if IsEof(line) then break; end if;
        keystr, valuestr := Splitter(line);
        data[KeyParser(keystr)] := FunctionOnLoad(ValueParser(valuestr));
    end while;
    return data;
    
end intrinsic;


intrinsic ReadCSV(fname : DataPath:=DatabaseDirectory(), FunctionOnLoad:=func<x|x>) -> Assoc
{}
    types := DetermineReadTypes(DataPath, fname);

    if #types eq 2 then
        return ReadCSVWithTypes(fname, types[1], types[2] :
                                DataPath:=DataPath, FunctionOnLoad:=FunctionOnLoad);
    end if;
    
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
            data[i] := FunctionOnLoad(args[2]);
        else
            data[i] := FunctionOnLoad(args[2..#args]);
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
    elif disc mod 8 eq 7 then
        return 1;
    else
        error "Unexpected discriminant residue.", disc mod 8;
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


intrinsic HalfWeil(counts::SeqEnum[FldRatElt]) -> SeqEnum[FldRatElt]
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

intrinsic IsSignAmbiguous(halfWeil) -> BoolElt, RngIntElt, RngIntElt
{From the first 11 coeffients of the Weil polynomial, 
determine if there is an ambiguity in the functional equation. Also return the
"Badness" rating -- the largest integer such that the coefficient of the degree 11+i term is 
zero.

If the sign is not ambiguous, return the sign as the third argument. The second argument
is set to -1 in this case (i.e., nonexistent badness).}

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
	sign := isWeil1 select 1 else -1;
        return false, -1, sign;
    end if;
end intrinsic;

intrinsic IsSignAmbiguous(list::SeqEnum[RngIntElt]) -> BoolElt, RngIntElt
{Given a list of point counts over F_2^k, k = 1,..., 11,
determine if there is an ambiguity in the functional equation. Also return the
"Badness" rating -- the largest integer such that the coefficient of the degree 11+i term is 
zero.

If the sign is not ambiguous, return the sign as the third argument. The second argument
is set to -1 in this case (i.e., nonexistent badness).}
    return IsSignAmbiguous(HalfWeil(list));    
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



intrinsic Charpoly(list::SeqEnum, sign::RngIntElt) -> RngUPolElt
{Given a list of point counts over F_2^k, k = 1,..., 11, and a sign for the functional equation,
returns the characteristic polynomial of Frobenius acting on primitive cohomology.}
    
    // p := t^22 + &+[c[i]*t^(22-i) : i in [1..11] ];
    // g := CONSTQt_ ! (t^22 * Evaluate(p, 1/t));

    c := HalfWeil(list);
    return WeilPolynomialFromHalf(c, sign);
end intrinsic;


intrinsic ZetaFunctionOfCubic(cubic :: RngMPolElt) -> RngUPolElt
{Given a smooth cubic form in six variables, return the primitive Weil polynomial.}
    counts := PointCounts(cubic);
    sign := FunctionalEquationSign(cubic);
    return WeilPolynomialFromHalf(HalfWeil(counts),  sign);
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


intrinsic AlgebraicRank(f :: RngUPolElt) -> RngUPolElt
{Given a Weil polynomia of the primitive cohomology of a cubic fourfold X, return the rank of the algebraic CH2(X).}
    return 23 - Degree(IrrationalFactor(f));
end intrinsic;




intrinsic GeometricRank(f :: RngUPolElt) -> RngUPolElt
{Given a Weil polynomial of the primitive cohomology of a cubic fourfold X, return the rank of the geometric CH2(X).}
    return 23- TranscendentalRank(f);
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
{Given the L polynomial of a cubic fourfold over F2, determine if the cubic fourfold has height 1 (i.e. is ordinary)}    
    g := Weil select CubicWeilPolynomialToLPolynomial(f) else f;
    NP := LowerVertices(NewtonPolygon(g, 2));

    hodgepolygon := [ <0, 0>, <1, 1>, <21, 41>, <22, 44>];
    return NP eq hodgepolygon;
end intrinsic;
