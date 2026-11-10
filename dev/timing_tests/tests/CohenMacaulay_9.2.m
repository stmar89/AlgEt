/*
    CohenMacaulay_9.2: CohenMacaulayType of all the 2288 overorders of the equation order E of f = (x^4+16)(x^4+81)
    (OverOrders not timed).
    Why: the order of WKICM_4.8 (Q(zeta_8) x Q(zeta_12), 6 singular primes with residue fields up to F_169): the
    types are at most 3, but there are 2288 trace dual ideals and indices to compute in an algebra of dimension 8.
    About 80 s in 2025, i.e. 35 ms per overorder: the per-order overhead of TraceDualIdeal/ColonIdeal/Index.
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=(x^4+16)*(x^4+81);
A:=EtaleAlgebra(f);
E:=EquationOrder(A);
oo:=OverOrders(E);
assert #oo eq 2288;
t0:=Cputime();
    cm:=[ CohenMacaulayType(T) : T in oo ];
t:=Cputime(t0);
assert Max(cm) eq 3;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
