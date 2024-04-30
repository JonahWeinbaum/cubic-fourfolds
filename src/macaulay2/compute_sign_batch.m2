load "Resultants.m2"

R = ZZ[u,v,w,x,y,z];
theseq = (); 
scanLines(l -> theseq = append(theseq, discriminant value l % 8), "input.m2");

print theseq;
exit;
-- "output.m2" << theseq << endl << close; --

