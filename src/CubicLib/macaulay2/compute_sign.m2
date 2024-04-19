load "Resultants.m2"

R = ZZ[u,v,w,x,y,z];
finput = value get "macaulay2/input.m2"

print discriminant(value(finput));
exit;
