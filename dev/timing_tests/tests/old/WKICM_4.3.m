/*
    WKICM_4.3: WKICM of the equation order E of f = x^8+16 (CM octic field, Weil polynomial for q = 2): 173 classes.
    Why: [O:E] = 2^14 and a unique singular prime, above 2, totally ramified in O; E is Gorenstein and WKICM uses the
    recursion over the minimal overorders; several overorders have Cohen-Macaulay type > 2. The 173 classes were
    first computed with the old method in 189000 s.
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=x^8+16;
A:=EtaleAlgebra(f);
R:=EquationOrder(A);
t0:=Cputime();
    wk:=WKICM(R);
t:=Cputime(t0);
assert #wk eq 173;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
