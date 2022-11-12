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
{ Given a sequence of 82 bytes, return the associated list of lines. The first 651 bits are 
interpreted as a boolean, where the i-th bit is 1 if the 'i'th line is in the associated list.
The remaining 5 zeroes for padding.

One further technicality -- we only guarantee 

    Set(lst) eq Set(DeserializeLinesThrough(SerializeLinesThrough(lst)))
}
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
