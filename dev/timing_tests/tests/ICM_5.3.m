/*
    ICM_5.3: ICM (GRH) of S = Z + 2 O in A = Q(sqrt(-5)) x Q(sqrt(-7)) x Q(sqrt(-11)), f = (x^2+5)(x^2+7)(x^2+11):
    3712 ideal classes.
    Why: three imaginary quadratic fields in which 2 behaves differently (ramified, split, inert), so O/2O has
    dimension 6 and the unique singular prime P = 2O of S has S/P = F_2; [O:S] = 2^5, S has Cohen-Macaulay type 5
    and 40 overorders, 28 of them not Gorenstein; Cl(O) = Z/2 and Pic(S) has 12 elements. WKICM(S) has 357 classes
    (the general branch, ~7 s) and the ideal class monoid 3712 (~13 s in 2025).
*/
AttachSpec(spec);
SetSeed(1);
SetClassGroupBounds("GRH");
_<x>:=PolynomialRing(Integers());
f:=(x^2+5)*(x^2+7)*(x^2+11);
A:=EtaleAlgebra(f);
O:=MaximalOrder(A);
S:=Order([One(A)] cat [2*z : z in ZBasis(O)]);
t0:=Cputime();
    icm:=ICM(S : GRH:=true);
t:=Cputime(t0);
assert #icm eq 3712;
assert #WKICM(S) eq 357;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
