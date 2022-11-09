
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
