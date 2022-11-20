//////////////////////////////////////////////
//
// LibPath.m
//
//////////////////////////////////////////////
//
// CHECK_TOKEN *Must* appear on Line 9. The EXACT formatting is important.
//
CHECK_TOKEN := "s;dljksgfpoierutqweoruimxj2385yu4q3ui5jrcmqp89742351n3c4958u1c";
//
// Utility to have a library discover its own path.
// The main usage is in case there are System calls
// to non-magma components that need to work independent
// of the directory.
// 

function ForceError()
    // Hmm... Using try-catch and e`Position, I can actually extract
    // the filename where the error occurs. I can likely hack this into
    // being able to get magma to set a path variable for system calls and
    // not worry about library location...Real skeevy, but should do the job.
    
    error "Forced Error from LibPath.m";
    return 1;
end function;

function PathFromErrorString(str)
    // Fish out the path from the error.
    S :=  Split(str, ",")[1];

    // Remove the file name from path.
    n := #"LibPath.m";
    filepath := S[10..#S-1]; 
    
    return filepath[1..#filepath-n];
end function;

// Standardize printing environment.
userCols := GetColumns();
userAutoCols := GetAutoColumns();
SetColumns(0);
SetAutoColumns(false);

// Extract Path.
try
    ForceError();
catch e
    CONST_LibPath := PathFromErrorString(e`Position);
    DEBUG_INFO := e;
end try;

// Reset User Column status.
SetColumns(userCols);
SetAutoColumns(userAutoCols);


// Verify that the path to this file has been correctly computed.
// If not, set the path value to the default empty string.
try
    // Extract Line 9 of THIS FILE.
    F := Open(CONST_LibPath * "LibPath.m", "r");
    for i in [1..8] do
	str := Gets(F);
    end for;
    
    // 16 is the length of "CHECK_TOKEN := \"";
    lineNine := Gets(F);
    checkStr := lineNine[17..#lineNine-2];
    error if checkStr ne CHECK_TOKEN, "Check Token comparison failure.";
catch e
    DEBUG_INFO_II := e;
    print "WARNING:", e`Position;
    print "Path to Libary cannot be automatically detected. Switching to default values.";
    print "Some functionality may not work.";
    CONST_LibPath := "";
end try;


//////////////////////////////
// Intrinsics

intrinsic PathToLib() -> MonStgElt
{Returns the directory of 'LibPath.m' as an absolute path.}
    return CONST_LibPath;
end intrinsic;

intrinsic PathToLibDebugInfo() -> MonStgElt
{Return debugging information.}
    if assigned DEBUG_INFO_II then
	return DEBUG_INFO, DEBUG_INFO_II, PathToLib();
    else
	return DEBUG_INFO, "Second value not assigned", PathToLib();
    end if;
end intrinsic;
