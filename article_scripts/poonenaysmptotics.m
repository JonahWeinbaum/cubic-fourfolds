//Computes the number of smooth hypersurfaces of degree d in P^(N-1) when N is not 


N := 3; //number of variables (dim X = N-2, X in P^(N-1))
d := 3; // d=degree(X)

//smoothfraction[d] = fraction of plane curves of degree d which are smooth

smoothfraction := AssociativeArray();

for d in [1..5] do
k := FiniteField(2);

R<[x]> := PolynomialRing(k, N);

V, bmap := GModule(GL(N,k), R, d);
P := ProjectiveSpace(R);
G := ActionGroup(V);

orbs := Orbits(G);
orbreps := [orbs[i][1] : i in [1..#orbs]];

//start at i=2 to toss out the zero orbit.
smoothorbsizes := [#orbs[i] : i in [2..#orbs] |IsNonsingular(Scheme(P, orbreps[i]@@bmap))];

smoothfraction[d] := &+smoothorbsizes/(#V-1);

end for;

//Results: smooth fraction, i.e. #U_d(F2)/(#P^(num. of monomials -1)(F2)). 


//Plane Curves(N=3)
//Poonen asymptotically says smoothfraction goes to 3/8 as d->infinity
//d=1: 1
//d=2: 4/9
//d=3: 112/341
//d=4: 1560/4681
//d=5: 98304/299593

//Surfaces in P3
//Poonen says smoothfraction goes to 21/64
//d=1: 1
//d=2: 448/1023
//d=3: 21504/69905
//d=4 (see quarticsurfaces.m): 10590854400/34359738367


//threefolds in P4
//Poonen says smoothfraction goes to 315/1024
//d=1: 1
//d=2:64/151
//d=3 (see bonus section of quarticsurfaces.m):330301440/1108378657 
//d=4: not feasible



//fourfolds in P5
//Poonen says smoothfraction goes to 9765/32768
//d=1: 1 
//d=2: 126976/299593
//d=3 (the paper): 1069562/ 3718649 

