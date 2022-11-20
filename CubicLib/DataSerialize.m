/////////////////////////////////////////////////////////
//
// Serialization and Deserialization functionality.
//
/////////////////////////////////////////////////////////

//Attach("DataRead.m");

intrinsic DeserializeLine(byteseq :: SeqEnum) -> ModMatFldElt
{'Given a sequence of 24 bits representing a line, return the associated 4 x 6 matrix. 
 The ordering of the entries is Magma\'s standard order.'}
        
    ret := [];
    seq := [];
    for byte in byteseq do
        seq cat:= [ BitwiseAnd(ShiftRight(byte, 8-i), 1): i in [1..8]];
    end for;

    return Matrix(GF(2), 4, 6, seq); 
end intrinsic;

intrinsic SerializeLine(form :: ModMatFldElt) -> SeqEnum
{'Given an echelon form, convert it into a sequence of bytes.'}
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


// Path to where linear subspace data is kept.
PATH_TO_LIB := PathToLib();
LINES_SUBSPACE_DIRECTORY := PATH_TO_LIB *
                            "../../database/linear_subspaces/lines_through_cubics/";
PLANES_SUBSPACE_DIRECTORY:= PATH_TO_LIB *
                            "../../database/linear_subspaces/planes_through_cubics/";


intrinsic ReadLinesIndex() -> SeqEnum
{Reads all lines in P5(F2) in a consistent ordering that will be used by calculations.
 Each line is stored as a matrix whose rows specify the linear equations cutting out the line.}
    lines := [];

    name := Sprintf(LINES_SUBSPACE_DIRECTORY * "lines.data");
    seq := [];
    file := Read(name);

    currloc := 1;

    while true do
       //Deserialize 3 Bytes and Advance File Pointer
    for i in [1..3] do
        byte := StringToCode(file[currloc]);
        seq cat:= [byte];
        currloc := currloc + 1;

    end for;
        Append(~lines, DeserializeLine(seq));
        seq := [];

    if currloc eq #file + 1 then
        break;
    end if;

    end while;
    return lines;
end intrinsic;




/*
intrinsic SerializeLinesThrough(linesthrough :: SeqEnum) -> SeqEnum
{Given a sequence of echelon forms representing lines through a cubic, convert it into a sequence of bytes.}
    byteseq := [];
    bitsseq := [];
    bits_num := Ceiling((#linesthrough + 1)*10 / 8) * 8;
    bytes_num := Ceiling((#linesthrough + 1)*10 / 8);

    //Add each line
    for i in [0..(#linesthrough - 1)] do
        for j in [1..10] do
            bitsseq[i*10 + j] := BitwiseAnd(ShiftRight(linesthrough[i+1], 10-j), 1);
        end for;
    end for;

    //Add delimiter (1024)
    for j in [1..10] do
        bitsseq[#(linesthrough)*10 + j] := 1;
    end for;
    //Add padding
    padding_start := (#linesthrough + 1)*10 + 1;
    for i in [padding_start..bits_num] do
        bitsseq[i] := 0;
    end for;

    //Convert to bytes sequence
    for i in [0..bytes_num-1] do 
        byte := 0;
        for j in [1..8] do
            current_bit := bitsseq[8*i + j];
            byte := BitwiseOr(byte, ShiftLeft(current_bit, 8-j));
end for;
        Append(~byteseq, byte);
    end for;

    return byteseq;
end intrinsic;

intrinsic DeserializeLinesThrough(byteseq :: SeqEnum) -> SeqEnum
{Given a sequence of bytes representing lines through a cubic, convert it into a sequence of indices.}
    linesthrough := [];
    n := 1;
    while Floor((n*10)/8) lt #(byteseq) do
        byte_num := Floor(n*10 / 8);
print(byte_num);
        line := 0;
        for i in [1..8] do 
            current_byte := byteseq[byte_num];
            current_bit := BitwiseAnd(ShiftRight(current_byte, 9-i), 1);
//print(current_bit);
            line := BitwiseOr(line, ShiftLeft(current_bit, 10-i));
        end for;
        for i in [1..2] do
            current_byte := byteseq[byte_num + 1];
            current_bit := BitwiseAnd(ShiftRight(current_byte, 9-i), 1);
//print(current_bit);
            line := BitwiseOr(line, ShiftLeft(current_bit, 10-i));
        end for;
//print(line);

        if line eq 1024 then 
            break;
        end if;

        Append(~linesthrough, line);
        
        n := n + 1;
    end while;

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
