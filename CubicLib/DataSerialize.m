Attach("CubicLib.m");
Attach("DataRead.m");
lines := ReadLinesIndex();

intrinsic SerializeLine(form :: ModMatFldElt) -> SeqEnum
{Given an echelon form, convert it into a sequence of bytes.}
    byteseq := [];
    for i in [0..2] do
        n := 1;
        bitstring := 0;
        for j in [1..8] do 
            if form[Ceiling(((j+1 + i*8) - 1) / 6)][(((j-1) + i*8) mod 6) + 1] eq 1 then
                bitstring := BitwiseOr(bitstring, ShiftLeft(1, 8-n));
            end if;
            n := n + 1;
        end for;
        Append(~byteseq, bitstring);
    end for;

    return byteseq;
end intrinsic;

intrinsic DeserializeLine(byteseq :: SeqEnum) -> ModMatFldElt
{Given a sequence of 24 bits representing a line, return the associated 4 x 6 matrix.
The ordering of the entries is Magma's standard order.}
        
    ret := [];
    seq := [];
    for byte in byteseq do
        seq cat:= [ BitwiseAnd(ShiftRight(byte, 8-i), 1): i in [1..8]];
    end for;

    return Matrix(GF(2), 4, 6, seq); 
end intrinsic;

intrinsic SerializeLinesThrough(linesthrough:: SeqEnum) -> SeqEnum
{Given a set of lines, convert it to a sequence of bytes.}
    byteseq := [];
    for i in [0..81] do
        bitstring := 0;
        bitstring := ShiftLeft(bitstring, 8);
        for j in [1..8] do 
            if i*8 + j in [652..656] then 
                bitstring := BitwiseOr(bitstring, ShiftLeft(0, 8-j));
                continue;
            end if;
            if lines[i*8 + j] in linesthrough then
                bitstring := BitwiseOr(bitstring, ShiftLeft(1, 8-j));
            end if;
        end for;
        Append(~byteseq, bitstring);
    end for;

    return byteseq;
end intrinsic;



intrinsic DeserializeLinesThrough(byteseq:: SeqEnum) -> SeqEnum
{Given a sequence of 82 bytes, with the first 651 bits each denoting whether the i'th line is
 contained in the cubic, and the remaining 5 zeroes for padding.}
    linesthrough := [];
    seq := [];
    for i in [0..81] do
        for j in [1..8] do
            if BitwiseAnd(byteseq[i+1], ShiftLeft(1, 8-j)) ne 0 then
                linesthrough cat:= [ lines[i*8 + j] ];
            end if;
        end for;
    end for;
    
    return linesthrough;
end intrinsic;
/*
intrinsic SerializePlane(form:: ModMatFldElt) -> SeqEnum
{Given an echeleon form, convert it to a sequence of bytes.}
    byteseq := [];
    for i in [0..2] do
        n := 1;
        bitstring := 0;
        if i eq 2 then 
            for j in [1..2] do 
                if form[Ceiling(((j+1 + i*8) - 1) / 6)][(((j-1) + i*8) mod 6) + 1] eq 1 then
                    bitstring := BitwiseOr(bitstring, ShiftLeft(1, 8-n));
                end if;
                n := n + 1;
            end for;
        else 
            for j in [1..8] do 
                if form[Ceiling(((j+1 + i*8) - 1) / 6)][(((j-1) + i*8) mod 6) + 1] eq 1 then
                    bitstring := BitwiseOr(bitstring, ShiftLeft(1, 8-n));
                end if;
                n := n + 1;
            end for;
        end if;
        Append(~byteseq, bitstring);
    end for;

    return byteseq;
end intrinsic;

intrinsic DeserializePlane(byteseq :: SeqEnum) -> ModMatFldElt
{Given a sequence of 18 bits representing a plane and 6 bits of zeroes for padding, return the
 associated 3 x6 matrix. The ordering of the entries is Magma's standard order.}
    ret := [];
    seq := [];
    n := 0;
    for byte in byteseq do
        if n eq 2 then 
            seq cat:= [ BitwiseAnd(ShiftRight(byte, 8-i), 1): i in [1..2]];
        else 
            seq cat:= [ BitwiseAnd(ShiftRight(byte, 8-i), 1): i in [1..8]];
        end if;
        n := n + 1;
    end for;
    
    return Matrix(k, 3, 6, seq); 
end intrinsic;

intrinsic SerializePlanesThrough(planesthrough:: SeqEnum) -> SeqEnum
{Given a set of planes, convert it to a sequence of bytes.}
    byteseq := [];
    for i in [0..174] do
        bitstring := 0;
        bitstring := ShiftLeft(bitstring, 8);
        for j in [1..8] do 
            if i*8 + j in [1396..1400] then 
                bitstring := BitwiseOr(bitstring, ShiftLeft(0, 8-j));
                continue;
            end if;
            if planes[i*8 + j] in planesthrough then
                bitstring := BitwiseOr(bitstring, ShiftLeft(1, 8-j));
            end if;
        end for;
        Append(~byteseq, bitstring);
    end for;

    return byteseq;
end intrinsic;

intrinsic DeserializePlanesThrough(byteseq:: SeqEnum) -> SeqEnum
{Given a sequence of 82 bytes, with the first 1395 bits each denoting whether the i'th plane is
 contained in the cubic, and the remaining 5 zeroes for padding.}
    planesthrough := [];
    seq := [];
    for i in [0..81] do
        for j in [1..8] do
            if BitwiseAnd(byteseq[i+1], ShiftLeft(1, 8-j)) ne 0 then
                planesthrough cat:= [ planes[i*8 + j] ];
            end if;
        end for;
    end for;
    
    return planesthrough;
end intrinsic;
*/
