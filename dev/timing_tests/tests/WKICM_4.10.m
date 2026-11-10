/*
    WKICM_4.10: WKICM of S = Z + 2 O in A = Q(zeta_8) x (the real cyclic cubic field of conductor 9),
    f = (x^4+1)(x^3-3x+1): 555 classes, 18 overorders.
    Why: a product of a CM field and a totally real field (signatures (0,2) and (3,0)); 2 is inert in the cubic field
    (f = 3) and totally ramified in Q(zeta_8) (e = 4), so the unique singular prime P = 2O has S/P = F_2 and
    (P:P)/P = O/2O = F_8 x F_2[t]/(t^4) of dimension 7; [O:S] = 2^6 and the Cohen-Macaulay type of S is 6 = n-1.
    Few overorders but many classes, found by the general branch of WKICM_bar through a large enumeration of subspaces
    over F_2 in dimension 7 (about 270 s in 2025).
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=(x^4+1)*(x^3-3*x+1);
A:=EtaleAlgebra(f);
O:=MaximalOrder(A);
S:=Order([One(A)] cat [2*z : z in ZBasis(O)]);
t0:=Cputime();
    wk:=WKICM(S);
t:=Cputime(t0);
assert #wk eq 555;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
