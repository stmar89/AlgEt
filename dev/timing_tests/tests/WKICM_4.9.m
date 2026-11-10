/*
    WKICM_4.9: WKICM of S = Z + 11 O in the cyclotomic field Q(zeta_5), f = x^4+x^3+x^2+x+1: 35 classes, 15 overorders.
    Why: 11 splits completely in Q(zeta_5), so the unique singular prime P = 11O has S/P = F_11 while
    (P:P)/P = O/11O = F_11^4; [O:S] = 11^3 and S has Cohen-Macaulay type 3 = n-1. The output is tiny, but the general
    branch of WKICM_bar enumerates subspaces of F_11-vector spaces and their orbits under (O/11O)^*: about 250 s in
    2025 for 35 classes, i.e. the cost of the enumeration with a large residue field, not of the output.
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=x^4+x^3+x^2+x+1;
A:=EtaleAlgebra(f);
O:=MaximalOrder(A);
S:=Order([One(A)] cat [11*z : z in ZBasis(O)]);
t0:=Cputime();
    wk:=WKICM(S);
t:=Cputime(t0);
assert #wk eq 35;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
