/*
    PicardGroup_1.3: PicardGroup (GRH) of the equation order E of f = (x^8+16)(x^8+81), a product of two CM octic
    fields (both factors are Weil polynomials, for q = 2 and q = 3).
    Why: a reducible f, so every step of PicardGroup works with direct sums (class groups, unit groups, ray residue
    rings); Cl(O) = Z/2. [O:E] ~ 2.8 10^24; E has 6 singular primes: above 2 (totally ramified in the first field),
    above 3 (e = 2 in the second), and two each above 5 and 13 with residue fields F_5^4 and F_13^4. This example
    used to spend ~8 s in ResidueField, which searched S/P for a root of a polynomial. Pic(E) has invariants
    [2,6,24,48,48,240].
*/
AttachSpec(spec);
SetSeed(1);
SetClassGroupBounds("GRH");
_<x>:=PolynomialRing(Integers());
f:=(x^8+16)*(x^8+81);
A:=EtaleAlgebra(f);
E:=EquationOrder(A);
t0:=Cputime();
    P,p:=PicardGroup(E : GRH:=true);
t:=Cputime(t0);
assert AbelianInvariants(P) eq [ 2, 6, 24, 48, 48, 240 ];
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
