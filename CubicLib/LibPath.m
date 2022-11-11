//////////////////////////////////////////////
//
// LibPath.m
//
//////////////////////////////////////////////
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

try
    ForceError();
catch e
    CONST_LibPath := PathFromErrorString(e`Position);
end try;

intrinsic PathToLib() -> MonStgElt
{Returns the path to the directory of 'LibPath.m'.}
    return CONST_LibPath;
end intrinsic;

