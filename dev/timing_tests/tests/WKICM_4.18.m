/*
    WKICM_4.18: WKICM of R = Z[F,q/F] for f = x^10-x^9+6x^8-10x^7+28x^6-48x^5+112x^4-160x^3+384x^2-256x+1024 (q = 4,
    degree 10: 44800 classes, 4608 overorders) and for f = x^8-x^7+10x^6-7x^5+138x^4-63x^3+810x^2-729x+6561 (q = 9,
    degree 8: 4880 classes, 1600 overorders).
    Why: the first f is (x^2-2x+4)(x^2+4)(x^2+3x+4)(x^4-2x^3+2x^2-8x+16), so A has four components (Q(sqrt(-3)), Q(i),
    Q(sqrt(-7)) and a quartic CM field), [O:R] = 2^9 3^11 5^2 and R has 4 singular primes; for the second f =
    (x^2-3x+9)(x^2+2x+9)(x^4-2x^2+81) (Q(sqrt(-3)), Q(sqrt(-2)) and a quartic CM field), [O:R] = 2^15 3 5^2 11^2 and R
    has 4 singular primes. WKICM splits the computation into the local orders R+P^kO, computes their classes by
    recursion and glues the local parts (the "divide and patch" branch). Found by a search over the Weil polynomials of
    degree <= 10 and q <= 9 for orders Z[F,q/F] of big index. Replaces WKICM_4.4 (Z[F,2/F] in degree 10 with 2 singular
    primes, 238 classes, about 3 s with 9362a5e).
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
// < f, q, number of weak equivalence classes of Z[F,q/F] >
data:=[ <x^10-x^9+6*x^8-10*x^7+28*x^6-48*x^5+112*x^4-160*x^3+384*x^2-256*x+1024,4,44800>,
        <x^8-x^7+10*x^6-7*x^5+138*x^4-63*x^3+810*x^2-729*x+6561,9,4880> ];
t:=0;
for d in data do
    f,q,n:=Explode(d);
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    R:=Order([F,q/F]);
    t0:=Cputime();
        wk:=WKICM(R);
    t+:=Cputime(t0);
    assert #wk eq n;
end for;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
