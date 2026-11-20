/*
    ICM_5.1: ICM (GRH) of the equation order E of f = (x^2+5)(x^2+20): 1356 ideal classes.
    Why: the two components are isomorphic to Q(sqrt(-5)) but given by different defining polynomials, so the package
    sees two distinct fields (Cl(O) = Z/2 x Z/2) that it cannot identify; [O:E] = 2 3^2 5^2 and E has 4 singular primes:
    two above 3 (3 splits in Q(sqrt(-5))), one above 2 and one above 5 (both ramified). E is Bass (quadratic
    components), so WKICM is just the 24 overorders; Pic(E) has 320 elements and the 1356 ideal classes come from
    the products of weak classes with Picard representatives, where the time goes (~12 s in 2025).
*/
AttachSpec(spec);
SetSeed(1);
SetClassGroupBounds("GRH");
_<x>:=PolynomialRing(Integers());
f:=(x^2+5)*(x^2+20);
A:=EtaleAlgebra(f);
E:=EquationOrder(A);
t0:=Cputime();
    icm:=ICM(E : GRH:=true);
t:=Cputime(t0);
assert #icm eq 1356;
assert #WKICM(E) eq 24 and #PicardGroup(E) eq 320;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
