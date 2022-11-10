


load "data-processing/read-data-planes-index.m";
load "data-processing/read-data-lines-index.m";
load "data-serialization/dataprocessing-planes.m";
load "data-serialization/dataprocessing-lines.m";

fourlines := [];

for i in  [1..#planes] do
l := [ln : ln in lines | RowSpace(planes[i]) subset RowSpace(ln) ];
fourlines[i] := {l[1], l[2], l[3], l[4]};
end for;

load  "data-processing/read-data-lines.m";

planesthrough := [];
time for i in [1..#linesthrough] do
    planesthrough[i] := [];
    for j in [1..1395] do
        if fourlines[j] subset linesthrough[i] then
            planesthrough[i] :=  planesthrough[i] cat [planes[j]];
        end if;
    end for;
    if (i mod 10000) eq 0 then print i;
    end if;
end for;
print "data computation finished";

fname := "../data/new/planes-data.data";
file := Open(fname, "w");

time for through in planesthrough do
WriteBytes(file, serialize_plane_info(through));
end for;
print "data file writing finished;"



/*
fourlineskeys := [ [Ar[ Matrix(w[i]) ] : i in [1..4]] : w in  fourlines];




planesdatakeyed := [];
for i in [1..85] do
time planesdatakeyed[i]:= [ [j : j in [1..1395] | (fourlineskeys[j][1] in ld)  and (fourlineskeys[j][2] in ld) and (fourlineskeys[j][3] in ld) and (fourlineskeys[j][4] in ld)   ]      : ld in linesdatakeyed[i]];
end for;*/
