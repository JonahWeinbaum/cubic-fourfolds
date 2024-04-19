//stabilizer_counts.csv actually records orbit size. This read-in script returns //the actual order of the stabilizers.

fname := "../database/group_action/stabilizers_info/stabilizer_counts.csv";
 stabsdata := AssociativeArray();

n := #GL(6, FiniteField(2));

 io := Open(fname, "r");
 print "stabilizer data loading...";
time while true do
 l := Gets(io);
 if IsEof(l) then break; end if;
 s := Split(l, ",");
 stabsdata[StringToInteger(s[1])] := n/StringToInteger(s[2]);
end while;

print "stabilizer data loaded."
