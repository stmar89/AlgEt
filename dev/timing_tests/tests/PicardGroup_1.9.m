/*
    PicardGroup_1.9: PicardGroup and UnitGroup (GRH) of the equation order E of f = (x^8+16)(x^8+81)(x^8+625), a product
    of three CM octic fields (the factors are Weil polynomials for q = 2, 3 and 5).
    Why: a reducible f with three components of dimension 8, so every step of PicardGroup works with direct sums (class
    groups, unit groups, residue rings) in dimension 24. The resultants of the factors make [O:E] = 2^54 3^20 5^20 7^8
    13^8 17^8 29^8 and E has 20 singular primes, so (O/f)^* is big: Pic(E) has 23 invariants, and the unit group of E
    has a huge index in the unit group of O. Replaces PicardGroup_1.3 (the same with (x^8+16)(x^8+81), about 1 s with
    9362a5e).
*/
AttachSpec(spec);
SetSeed(1);
SetClassGroupBounds("GRH");
_<x>:=PolynomialRing(Integers());
f:=(x^8+16)*(x^8+81)*(x^8+625);
A:=EtaleAlgebra(f);
assert [ DefiningPolynomial(K) : K in Components(A) ] eq [ x^8+16, x^8+81, x^8+625 ];
E:=EquationOrder(A);
t0:=Cputime();
    P,p:=PicardGroup(E : GRH:=true);
    U,u:=UnitGroup(E : GRH:=true);
t:=Cputime(t0);
assert AbelianInvariants(P) eq [ 2, 2, 4, 4, 8, 8, 16, 16, 16, 48, 48, 48, 48, 48, 48, 48, 48, 480, 480, 480, 480,
                                 10080, 10080 ];
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
