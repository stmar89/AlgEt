/*
    WKICM_4.4: WKICM of R = Z[F,2/F] for f = x^10-x^9+4x^8-6x^7+8x^6-16x^5+16x^4-24x^3+32x^2-16x+32, a Weil
    polynomial for q = 2 of degree 10: 238 classes.
    Why: f = (x^2-x+2)(x^2+2)(x^2+2x+2)(x^4-2x^3+2x^2-4x+4): the algebra has 4 components (Q(sqrt(-7)), Q(sqrt(-2)),
    Q(i) and a quartic CM field), [O:R] = 93312 = 2^7 3^6, and R has 2 singular primes: one above 3 with residue
    field F_9 and one above 2 with residue field F_2. WKICM splits the computation into the two local orders R+P^kO,
    computes their classes by recursion, and glues the local parts (the "divide and patch" branch).
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=x^10 - x^9 + 4*x^8 - 6*x^7 + 8*x^6 - 16*x^5 + 16*x^4 - 24*x^3 + 32*x^2 - 16*x + 32;
A:=EtaleAlgebra(f);
F:=PrimitiveElement(A);
q:=Integers() ! Round(ConstantCoefficient(f)^(2/Degree(f)));
R:=Order([F,q/F]);
t0:=Cputime();
    wk:=WKICM(R);
t:=Cputime(t0);
assert #wk eq 238;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
