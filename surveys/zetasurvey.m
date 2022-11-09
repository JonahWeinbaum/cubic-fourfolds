// Host/port data.
host := "localhost";
port := 10000;
Attach("CubicLib.m");


AttachSpec("~/magma-parallel-cookbook/spec");
socket := Socket(: LocalHost := "localhost", LocalPort := 10000);

// Load database
orbdata := LoadCubicOrbitData(: Flat:=true, Quick:=false); // 2 minute load.

// Launch!
orbzip := [<orbdata[i], i> : i in [1..#orbdata]];

StartDistributedWorkers("zeta-worker.m", 40);
results := DistributedManager(socket, orbzip);

print "";

// Resource cleanup.
delete socket;

