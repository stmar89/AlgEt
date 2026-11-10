/*
    OverOrders_3.8: OverOrders of S = Z + 2 O in the cyclotomic field Q(zeta_16), f = x^8+1: 110 overorders.
    Why: O = Z[zeta_16], 2 is totally ramified (e = 8) and [O:S] = 2^7; S has a unique singular prime P = 2O with
    S/P = F_2 and (P:P)/P = O/2O of dimension 7, and Cohen-Macaulay type 7 = n-1, the maximum possible. The local
    lattice of overorders at P is deep: 74 of the 110 overorders are not Gorenstein. A test of MinimalOverOrdersAtPrime
    with a large T/P (subspaces of F_2^7 and their Frobenius-invariance), where the type-2 enumeration is the cost.
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=x^8+1;
A:=EtaleAlgebra(f);
O:=MaximalOrder(A);
S:=Order([One(A)] cat [2*z : z in ZBasis(O)]);
t0:=Cputime();
    oo:=OverOrders(S);
t:=Cputime(t0);
assert #oo eq 110;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
