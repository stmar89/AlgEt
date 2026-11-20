/*
    WKICM_4.6: WKICM of R = Z[F,25/F] for f = x^6+8x^5+50x^4+200x^3+1250x^2+5000x+15625, a Weil polynomial for
    q = 25 (three imaginary quadratic fields, see OverOrders_3.3): 560 classes.
    Why: [O:R] = 11407500 = 2^2 3^2 5^5 13^2; R has 5 singular primes (above 2 with residue field F_4, 3, 5 and two
    above 13), the most among the examples, so the gluing step of WKICM dominates: the local parts are small, but 560
    global classes have to be patched together from 5-tuples of local classes, each with intersections of prime
    powers and LLL reductions. R is Gorenstein; compare with WKICM_4.2 (an overorder of type 3 in the same algebra)
    and OverOrders_3.3 (the equation order).
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=x^6 + 8*x^5 + 50*x^4 + 200*x^3 + 1250*x^2 + 5000*x + 15625;
A:=EtaleAlgebra(f);
F:=PrimitiveElement(A);
R:=Order([F,25/F]);
t0:=Cputime();
    wk:=WKICM(R);
t:=Cputime(t0);
assert #wk eq 560;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
