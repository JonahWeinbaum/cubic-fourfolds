// compute_zeta_functions.m
//
// To compute all the zeta functions efficiently, one should *not* use the folloing code:
//
//    for f in SMOOTH_CUBIC_REPRESENTATIVES do
//        ZetaFunctionOfCubic(f);
//        <Write to file>
//    end for;
//
// The reason is that while this is convenient for one-off calls, the overhead of launching
// sage or macaulay2 is prohibitive. Computing the result for all ~1.06 million cubics would
// take far longer than it should. This script batches the calls to Macaulay2 and is much
// faster.
AttachSpec("../CubicLib/CubicLib.spec");

function FormatForM2(f)
    R<x,y,z,u,v,w> := PolynomialRing(BaseRing(Parent(f)), 6);
    fR := R ! f;
    fstr := Sprint(fR);
    return fstr;
end function;

function ParseM2List(s)
    if #s eq 0 then error "Empty output recieved from M2."; end if;

    seq := false;
    for S in Split(s, "()") do
	if Regexp("^[0-9 ,]*$", S) then
	    seq := S;
	    break;
	end if;
    end for;

    if seq cmpeq false then
	error "Parsing error on input: ", s;
    end if;

    list := "[" * seq * "]"; // Replace parens.
    return eval list;
end function;

function FESign(disc)
    // Once the discriminant is found, compute the sign.
    if disc mod 8 eq 3 then
        return -1;
    elif disc mod 8 eq 7 then
        return 1;
    else
        error "Unexpected discriminant residue.", disc mod 8;
    end if;
end function;


////////////////////////////////////////
//
// Main script.
//
////////////////////////////////////////

// 1. Load data.
OUTPUT_FILE := DatabaseDirectory() * "zeta_functions/zeta_coefficients.csv";

time cubics := LoadCubicOrbitData(: Flat);
time Counts := ReadCSV("zeta_functions/point_counts.csv");
time smoothIndices := Setseq(Keys(ReadCSV("differentiability/smooth/smooth.csv")));

// 2. Now figure out a good enough batch size. Our Macaulay2 script is not memory efficient,
//    so one must be careful not to make the batches too large.

BATCH_SIZE := 1000;
NUM_CUBICS := #smoothIndices;

istart := 1;
while istart le NUM_CUBICS do
    iend := istart + BATCH_SIZE - 1;
    iend := (iend gt NUM_CUBICS) select NUM_CUBICS else iend;

    // Write the batch to file.
    outputString := "";
    for j in smoothIndices[istart..iend-1] do
	f := cubics[j];
	outputString cat:= FormatForM2(f) * "\n";
    end for;
    batchlast := smoothIndices[iend];
    outputString cat:= FormatForM2(cubics[batchlast]); // No ending newline.
    Write("input.m2", outputString : Overwrite);

    // Run M2.
    cmd := Sprint("M2 --script ../macaulay2/compute_sign_batch.m2");
    time m2output := Read(POpen(cmd, "r"));
    discs := ParseM2List(m2output);

    // Write zeta functions to file..
    for j->cindex in smoothIndices[istart..iend] do
	Lf := Charpoly(Counts[cindex], FESign(discs[j]));
	line := Sprintf("%o, %o", cindex, Coefficients(Lf));
	Write(OUTPUT_FILE, line);
    end for;
    
    // Increment.
    istart := iend + 1;
end while;

// Cleanup temp files.
dont_print := Read(POpen("rm input.m2", "r"));
