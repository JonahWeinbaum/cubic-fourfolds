R.<u,v,w,x,y,z> = PolynomialRing(ZZ, 6);

f = R(f)
D = R.macaulay_resultant([f.derivative(m) for m in R.gens()])

print(D)
quit;

