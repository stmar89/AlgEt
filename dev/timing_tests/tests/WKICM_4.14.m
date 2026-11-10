/*
    WKICM_4.14: WKICM of the equation order E of f = x^10-x^9+4x^8-6x^7+8x^6-16x^5+16x^4-24x^3+32x^2-16x+32, a Weil
    polynomial for q = 2 of degree 10: 12096 classes, 2286 overorders.
    Why: 4 components (Q(sqrt(-7)), Q(sqrt(-2)), Q(i) and a quartic CM field), [O:E] = 2^17 3^6, 2 singular primes
    (above 3 with residue field F_9, above 2 with four primes of O above it); E is Gorenstein, Pic(E) has 3072 elements.
    The largest number of weak classes in the suite after WKICM_4.7 (about 320 s in 2025): a stress test of the
    recursion over 2286 overorders, of the gluing of two local parts, and of the Hash-based sets of ideals. Compare
    with GraphOverOrders_3.9 (the same order).
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=x^10 - x^9 + 4*x^8 - 6*x^7 + 8*x^6 - 16*x^5 + 16*x^4 - 24*x^3 + 32*x^2 - 16*x + 32;
A:=EtaleAlgebra(f);
E:=EquationOrder(A);
t0:=Cputime();
    wk:=WKICM(E);
t:=Cputime(t0);
assert #wk eq 12096;
assert #OverOrders(E) eq 2286;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
