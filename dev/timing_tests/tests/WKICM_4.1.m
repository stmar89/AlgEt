/*
    WKICM_4.1: WKICM of R = Z[F,2/F] in the CM octic field of f = x^8+16 (F the primitive element), after
    OverOrders(R : populateoo_in_oo) (not timed): 31 weak equivalence classes, 19 overorders; computed on 10 fresh
    copies of the algebra (the result is cached in the order).
    Why: f is a Weil polynomial for q = 2 (supersingular), R is the order Z[F,q/F] that occurs for abelian varieties
    (it contains the Verschiebung q/F). [O:R] = 256; the unique singular prime is above 2, which is totally ramified
    in O (e = 8); R is Gorenstein but not Bass (some overorders have Cohen-Macaulay type > 1). A quick test of the
    WKICM recursion over the overorders.
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=x^8+16;
t:=0;
for i in [1..10] do
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    R:=Order([F,2/F]);
    oo:=OverOrders(R : populateoo_in_oo:=true);
    assert #oo eq 19;
    t0:=Cputime();
        wk:=WKICM(R);
    t+:=Cputime(t0);
    assert #wk eq 31;
end for;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
