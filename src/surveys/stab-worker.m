// Same as the manager process.
host := "localhost";
port := 10000;
AttachSpec("../CubicLib/CubicLib.spec");

// load "computecharpoly.m"; // 2 minutes load.

function ReportOrbitSizes(tuple)

    index := tuple[2];
    f := tuple[1];
    
    try
	a := CubicOrbitSize(f);
	ok_write := WriteOrbitSizeData(index, a);
    catch e
	ReportError(index, e);
    end try;
    return 0;
end function;

// Activate the worker
DistributedWorker(host, port, ReportOrbitSizes);

// Terminate the worker once it is done.
quit;
