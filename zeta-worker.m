// Same as the manager process.
host := "localhost";
port := 10000;
Attach("CubicLib.m");

load "computecharpoly.m"; // 2 minutes load.


function ReportPointCounts(cubic_blist)

    f := b @@ Bit;
    try
	a := PointCounts(f);
	WriteZetaData(i, 1, IsSmooth(f), a);
    catch e
	print e;
    end try;
    return 0;
end function;

// Activate the worker
DistributedWorker(host, port, ReportPointCounts);

// Terminate the worker once it is done.
quit;
