/*
    WKICM_4.7: WKICM of the equation order of f = x^8-2x^7+7x^6-14x^5+40x^4-56x^3+112x^2-128x+256, a Weil polynomial
    for q = 2 of degree 8: 114492 classes. About 22 minutes on a laptop (2025): the biggest example of the package.
    Why: a single singular prime with a huge local lattice of overorders and many overorders of Cohen-Macaulay type
    > 2; this is the stress test for the recursion of WKICM, for the Hash-based sets of orders and ideals, and for the
    memory footprint. v1.0.2 does not finish it within the default 40-minute limit (aerox, babbage, laplace).
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=x^8 - 2*x^7 + 7*x^6 - 14*x^5 + 40*x^4 - 56*x^3 + 112*x^2 - 128*x + 256;
A:=EtaleAlgebra(f);
R:=EquationOrder(A);
t0:=Cputime();
    wk:=WKICM(R);
t:=Cputime(t0);
assert #wk eq 114492;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
