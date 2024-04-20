// Host/port data.
host := "localhost";
port := 10000;

AttachSpec("../CubicLib/CubicLib.spec");
AttachSpec("~/magma-parallel-cookbook/spec");
socket := Socket(: LocalHost := "localhost", LocalPort := 10000);

// Load database
orbdata := LoadCubicOrbitData(: Flat:=true, Quick:=false); // 2 minute load.

// Launch!
orbzip := [<orbdata[i], i> : i in [1..#orbdata]];

StartDistributedWorkers("stab-worker.m", 40);
results := DistributedManager(socket, orbzip);

print "";

// Resource cleanup.
delete socket;

