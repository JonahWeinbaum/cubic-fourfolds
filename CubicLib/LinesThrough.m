/////////////////////////////////////////////////
//
// LinesThrough.m
//
/////////////////////////////////////////////////

// Technically, this is why it needs to be in a separate file. Otherwise Magma doesn't
// like assigning this constant.
CONST_lines := ReadLinesIndex();

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
            if CONST_lines[i*8 + j] in linesthrough then
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
	    linesthrough cat:= [ CONST_lines[i*8 + j] ];
            end if;
        end for;
    end for;
    
    return linesthrough;
end intrinsic;
