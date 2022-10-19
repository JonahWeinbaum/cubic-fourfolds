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

DATA_DIRECTORY := "../data/zeta/";
CUBIC_ID_FILE := "orbrep.csv";
ISSMOOTH_FILE := "smooth.csv";
POINT_COUNTS_FILE := "point_counts.csv";
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

intrinsic LoadCubicOrbitData(: RemoveZero:=true, Flat:=false, Quick:=false) -> SeqEnum
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
	     name2 := Sprintf("../data/new/orbitreps-%o.data", k);
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
    
intrinsic Charpoly(list::SeqEnum[RngIntElt]) -> RngUPolElt
    {input list of point counts over F_2^k, k = 1,..., 11, returns charpoly of frobenius on nonprimitive cohomology.}
    tr := [(list[m] - 1 - 2^m - 4^m - 8^m - 16^m)/4^m : m in [1..11]];
    c := [];
    c[1] := -tr[1];
    for k in [2..11] do
        c[k] := (tr[k] + &+[c[i]*tr[k-i] : i in [1..k-1]] )/(-k);
    end for;

    p := t^22 + &+[c[i]*t^(22-i) : i in [1..11] ];
    g := CONSTQt_ ! (t^22 * Evaluate(p, 1/t));

       // TODO: fix this.
    if c[11] ne 0 then
        ret := p + Modexp(g, 1, t^11);
        return ret;
    else
        poly1 := p + Modexp(g, 1, t^11);
        poly2 := p - Modexp(g, 1, t^11);
        roots1 := Roots(Evaluate(poly1, s));
        roots2 := Roots(Evaluate(poly2, s));
        mods1 := [Modulus(roots1[i][1]) : i in [1..#roots1]];
        mods2 := [Modulus(roots2[i][1]) : i in [1..#roots2]];

    one := ComplexField(30)! 1;
    if [Integers()! mods1[i] eq one : i in [1..#mods1] ] eq [true  : i in [1..#mods1] ] then
        return poly1;
        else return poly2;
    end if;

    end if;
    
end intrinsic;
        
