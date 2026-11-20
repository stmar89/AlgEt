/*
    WKICM_4.22: WKICM of S = Z + 4O in the totally real sextic field of f = x^6-6x^4+9x^2-3: 3382 classes, 502
    overorders.
    Why: the totally real example of the group (signature (6,0): the unit group has rank 5); 2 has e = 2 and f = 3 in
    O, so O/4O has 2^12 elements; the unique singular prime P = 2Z + 4O has S/P = F_2, [O:S] = 2^10 and S has
    Cohen-Macaulay type 5 = n-1. Compared with Z + 2O the local lattice is two layers deeper, and WKICM_bar runs on
    hundreds of non-Gorenstein overorders. Replaces WKICM_4.13 (Z + 2O in the same field, 66 classes, about 3 s with
    9362a5e).
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=x^6-6*x^4+9*x^2-3;
A:=EtaleAlgebra(f);
O:=MaximalOrder(A);
S:=Order([One(A)] cat [4*z : z in ZBasis(O)]);
t0:=Cputime();
    wk:=WKICM(S);
t:=Cputime(t0);
assert #wk eq 3382;
assert #OverOrders(S) eq 502;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
