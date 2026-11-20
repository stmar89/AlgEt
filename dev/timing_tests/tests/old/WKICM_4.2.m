/*
    WKICM_4.2: WKICM_bar of an order S, given by a Z-basis, in the algebra of f = x^6+8x^5+50x^4+200x^3+1250x^2+
    5000x+15625 (three imaginary quadratic fields, see OverOrders_3.3): 7 weak equivalence classes with multiplicator
    ring S, computed on 6 fresh copies of the algebra (the result is cached in the order).
    Why: [O:S] = 2281500 = 2^2 3^2 5^3 13^2; S has 5 singular primes (above 2 with residue field F_4, 3, 5 and two
    above 13) and Cohen-Macaulay type 3, so WKICM_bar goes through the general branch: local computations at each
    singular prime (orders S+P^kO), then the gluing of the local classes.
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=x^6 + 8*x^5 + 50*x^4 + 200*x^3 + 1250*x^2 + 5000*x + 15625;
dataS:=[
<[1,0],[1,0],[1,0]>,
<[0,1],[0,1],[0,1]>,
<[0,-2],[0,4],[0,-2]>,
<[8,0],[-2,0],[-5,0]>,
<[0,-6],[0,-2],[0,7]>,
<[2,0],[-18,0],[15,0]>
];
t:=0;
for i in [1..6] do
    A:=EtaleAlgebra(f);
    S:=Order([ A ! g : g in dataS ]);
    t0:=Cputime();
        wk:=WKICM_bar(S);
    t+:=Cputime(t0);
    assert #wk eq 7;
end for;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
