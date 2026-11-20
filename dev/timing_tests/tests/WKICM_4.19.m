/*
    WKICM_4.19: WKICM of R = Z[F,7/F] for f = x^10-4x^8-5x^6-35x^4-1372x^2+16807, a Weil polynomial for q = 7 of degree
    10: 28288 classes, 6860 overorders.
    Why: f = (x^2+7)(x^4-13x^2+49)(x^4+2x^2+49): Q(sqrt(-7)) and two quartic CM fields; f is not ordinary. [O:R] = 2^11
    3^14 5^4 and R has 4 singular primes, with deep local lattices above 2 and 3: the classes are the product of the
    sets of local classes. Found by a search over the Weil polynomials of degree <= 10 and q <= 9 for orders Z[F,q/F] of
    big index. Replaces WKICM_4.5 (Z[F,25/F] for a sextic, 3 singular primes, 315 classes, about 2 s with 9362a5e).
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=x^10-4*x^8-5*x^6-35*x^4-1372*x^2+16807;
q:=7;
A:=EtaleAlgebra(f);
F:=PrimitiveElement(A);
R:=Order([F,q/F]);
t0:=Cputime();
    wk:=WKICM(R);
t:=Cputime(t0);
assert #wk eq 28288;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
