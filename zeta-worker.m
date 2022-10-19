// Same as the manager process.
host := "localhost";
port := 10000;
Attach("CubicLib.m");

load "computecharpoly.m"; // 2 minutes load.


function ReportPointCounts(tuple)

    index := tuple[2];
    blist := tuple[1];
    
    f := BitListToCubic(blist);
    try
	a := PointCounts(f, index);
	WriteZetaData(index, IsSmooth(f), a);
    catch e
	ReportError(index, e);
    end try;
    return 0;
end function;

// Activate the worker
DistributedWorker(host, port, ReportPointCounts);

// Terminate the worker once it is done.
quit;
