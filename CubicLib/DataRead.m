Attach("CubicLib.m");
// Attach("DataProcessing.m");

// Path to where linear subspace data is kept.
LINES_SUBSPACE_DIRECTORY := "../../database/linear_subspaces/lines_through_cubics/";
PLANES_SUBSPACE_DIRECTORY := "../../database/linear_subspaces/planes_through_cubics/";

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
intrinsic ReadLines() -> SeqEnum
{Reads data of lines through each cubic. Each line is stored as a matrix whose rows specify 
the linear equations cutting out the line.}
    lines := [];

    name := Sprintf(LINES_SUBSPACE_DIRECTORY * "lines-through-indexed.data");
    linesthrough := [];
    file := Read(name);

    currloc := 1;
    n := 1;
    while true do

	//Deserialize 3 Bytes and Advance File Pointer 
	for i in [1..82] do
	    byte := StringToCode(file[currloc]);
	    seq cat:= [byte];
	    currloc := currloc + 1;
            
	end for;
	Append(~linesthrough, DeserializeLinesThrough(seq));
	seq := [];


	if currloc eq #file + 1 then
	    break;
	end if; 
    end while;
    return lines;
end intrinsic;

intrinsic ReadPlanesIndex() -> SeqEnum
{Reads all planes in P5(F2) in a consistent ordering that will be used by calculations. 
Each plane is stored as a matrix whose rows specify the linear equations cutting out the plane.}

    planes := [];

    name := Sprintf(PLANES_SUBSPACE_DIRECTORY * "planes-index.data");
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
	Append(~planes, DeserializePlane(seq));
	seq := [];

	if currloc eq #file + 1 then
	    break;
	end if; 
    end while;
    return planes;
end intrinsic;
*/
