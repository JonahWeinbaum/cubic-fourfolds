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
	is_smth := IsSmooth(f);
	if is_smth then
	    a := PointCounts(f : ExecNum := index);
	    ok_write := WriteZetaData(index, is_smth, a);
	end if;
    catch e
	ReportError(index, e);
    end try;
    return 0;
end function;

// Activate the worker
DistributedWorker(host, port, ReportPointCounts);

// Terminate the worker once it is done.
quit;
