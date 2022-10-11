// Host/port data.
host := "localhost";
port := 10000;
Attach("CubicLib.m");

// Set test flag.
PARALLEL_MODE := false;

if PARALLEL_MODE then
    AttachSpec("~/magma-parallel-cookbook/spec");
    socket := Socket(: LocalHost := "localhost", LocalPort := 10000);

    // Load database
    orbdata := LoadCubicOrbitData(: Flat:=true); // 2 minute load.

    // Launch!
    StartDistributedWorkers("zeta-worker.m", 40);
    results := DistributedManager(socket, orblist);

    print "";
    
    // Resource cleanup.
    delete socket;
    
else
    // Load Point counting function.
    // load "computecharpoly.m";
    orbdata := LoadCubicOrbitData(: Flat:=false); // 2 minute load.

    for i in [1..85] do
        for j in [1..1] do
            b := orbdata[i][j]; // For example;

            try
	        f := BitListToCubic(b);
	        a := PointCounts(f);
	        print a; // Want to parallelize this call.
	        WriteZetaData(i, j, IsSmooth(f), a);
            catch e
	        print e;
            end try;
        end for;
    end for;
    
end if;

