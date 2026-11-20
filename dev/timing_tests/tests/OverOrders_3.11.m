/*
    OverOrders_3.11: OverOrders of the equation order E = Z[2 zeta_15] of f = 2^8 Phi_15(x/2) = x^8-2x^7+8x^5-16x^4+
    32x^3-128x+256, in the cyclotomic field Q(zeta_15): 6904 overorders.
    Why: [O:E] = 2^28 and E has a unique singular prime P, above 2, with E/P = F_2, while 2 is unramified in O and
    splits into two primes of residue degree 4. A single deep local lattice: the time goes into
    MinimalOverOrdersAtPrime on thousands of orders. Replaces OverOrders_3.2 (3 singular primes, 3312 overorders,
    about 2 s with 9362a5e).
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=x^8-2*x^7+8*x^5-16*x^4+32*x^3-128*x+256;
A:=EtaleAlgebra(f);
E:=EquationOrder(A);
t0:=Cputime();
    oo:=OverOrders(E);
t:=Cputime(t0);
assert #oo eq 6904;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
