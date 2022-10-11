// Host/port data.
host := "localhost";
port := 10000;
Attach("CubicLib.m");

AttachSpec("~/magma-parallel-cookbook/spec");
socket := Socket(: LocalHost := "localhost", LocalPort := 10000);



// Load database.
//load "read-data-eof.m"; 
orbdata := LoadCubicOrbitData(: Flat:=true); // 2 minute load.

// Load the bit conversion function;

// Loads orbdata, an array of 85 arrays.
// Each entry of the inner array is a length 56 bit string, 

//orblist := [f : f in &cat orbdata | Seqint(f,2) ne 0];

//StartDistributedWorkers("zeta-worker.m", 40);
//results := DistributedManager(socket, orblist);


// Print something to get the prompt back.
print "";

// b @@ bit to converts to cubic.

// Check if it is smooth.

// Loop PointsCounts over all smooth orbit classes.


for i in [1..85] do
    b := orbdata[i]; // For example;

    try
	f := BitListToCubic(b);
	a := PointCounts(f);
	print a; // Want to parallelize this call.
	WriteZetaData(i, 1, IsSmooth(f), a);
    catch e
	print e;
    end try;
end for;



// Resource cleanup.
delete socket;
