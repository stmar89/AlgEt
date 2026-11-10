/*
    WKICM_4.11: WKICM of the equation order E of A = K x K with K = Q(zeta_8) given by x^4+16: 3468 classes, 1452
    overorders.
    Why: the two components are the same number field object (DefiningPolynomial(A) is not defined); Cl(O) is trivial,
    [O:E] = 2^12 3^4 7^2 257 and E has 5 singular primes: two above 2 (one per copy, totally ramified), one each above 3
    and 7 with residue fields F_9 and F_49, and one above 257, which splits. E is Gorenstein and Pic(E) has 12288
    elements. Compare with WKICM_4.8 (two different quartic CM fields, similar size).
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
K:=NumberField(x^4+16);
A:=EtaleAlgebra([K,K]);
E:=EquationOrder(A);
t0:=Cputime();
    wk:=WKICM(E);
t:=Cputime(t0);
assert #wk eq 3468;
assert #OverOrders(E) eq 1452;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
