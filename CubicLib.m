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

intrinsic LoadCubicOrbitData(: RemoveZero:=true, Flat:=false) -> SeqEnum
{Loads the precomputed orbit data.}

    ZEROCUBIC := [0 : i in [1..56]];
    print "loading data...";

    orbdata := [];
    time for k in [1..85] do
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

intrinsic WriteZetaData(i, j, issmooth, pointcounts) -> RngIntElt
{}
    Write(DATA_DIRECTORY * CUBIC_ID_FILE, Sprintf("%o,%o", i, j));    
    Write(DATA_DIRECTORY * ISSMOOTH_FILE, issmooth);
    Write(DATA_DIRECTORY * POINT_COUNTS_FILE, pointcounts);    
    return 0;
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

