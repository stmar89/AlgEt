/*
    WKICM_4.13: WKICM of S = Z + 2 O in the totally real sextic field of f = x^6-6x^4+9x^2-3: 66 classes, 18 overorders.
    Why: the only totally real example of the group (signature (6,0): the unit group has rank 5, the class group of O
    is trivial, Pic(S) is trivial); 2 has e = 2 and f = 3 in O, so the unique singular prime P = 2O has S/P = F_2 and
    (P:P)/P = O/2O of dimension 5; [O:S] = 2^5 and S has Cohen-Macaulay type 5 = n-1, with 9 of the 18 overorders not
    Gorenstein. The general branch of WKICM_bar in a totally real setting, in ~3 s.
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=x^6-6*x^4+9*x^2-3;
A:=EtaleAlgebra(f);
O:=MaximalOrder(A);
S:=Order([One(A)] cat [2*z : z in ZBasis(O)]);
t0:=Cputime();
    wk:=WKICM(S);
t:=Cputime(t0);
assert #wk eq 66;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
