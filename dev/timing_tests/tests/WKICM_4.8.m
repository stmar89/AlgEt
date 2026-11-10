/*
    WKICM_4.8: WKICM of the equation order E of f = (x^4+16)(x^4+81): 4896 classes, 2288 overorders.
    Why: two CM quartic fields, Q(zeta_8) and Q(zeta_12) (the factors are Weil polynomials for q = 2 and q = 3), with
    trivial class groups; [O:E] = 2^6 3^6 5^4 13^4 and E has 6 singular primes: one above 2 (totally ramified in the
    first field), one above 3, and two each above 5 and 13 with residue fields F_25 and F_169. E is Gorenstein and
    Pic(E) has 663552 elements. Many small local lattices and a big product: the gluing of the local weak classes
    ("divide and patch") dominates.
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=(x^4+16)*(x^4+81);
A:=EtaleAlgebra(f);
E:=EquationOrder(A);
t0:=Cputime();
    wk:=WKICM(E);
t:=Cputime(t0);
assert #wk eq 4896;
assert #OverOrders(E) eq 2288;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
