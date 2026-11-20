/*
    WKICM_4.20: WKICM of R = Z[F,5/F] for f = x^10+2x^9+2x^8-8x^7-11x^6-36x^5-55x^4-200x^3+250x^2+1250x+3125, an
    ordinary Weil polynomial for q = 5 of degree 10: 68736 classes, 20352 overorders.
    Why: f = (x^2-2x+5)(x^2+3x+5)(x^2+4x+5)(x^4-3x^3+4x^2-15x+25): Q(i) twice, given by different polynomials,
    Q(sqrt(-11)) and a quartic CM field; [O:R] = 2^17 3^4 5^2 11^2 and R has 8 singular primes, so the gluing step of
    WKICM dominates: the local parts are small, but 68736 global classes are patched together from 8-tuples of local
    classes, each with intersections of prime powers and LLL reductions. Found by a search over the Weil polynomials of
    degree <= 10 and q <= 9 for orders Z[F,q/F] of big index. Replaces WKICM_4.6 (Z[F,25/F] in three imaginary quadratic
    fields, 5 singular primes, 560 classes, about 3 s with 9362a5e).
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=x^10+2*x^9+2*x^8-8*x^7-11*x^6-36*x^5-55*x^4-200*x^3+250*x^2+1250*x+3125;
q:=5;
A:=EtaleAlgebra(f);
F:=PrimitiveElement(A);
R:=Order([F,q/F]);
t0:=Cputime();
    wk:=WKICM(R);
t:=Cputime(t0);
assert #wk eq 68736;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
