/*
    OverOrders_3.7: OverOrders of the equation order E of A = K x K, where K is the CM quartic field defined by the
    Weil polynomial x^4+x^2+529 (q = 23): 36608 overorders.
    Why: the two components are the same number field object (DefiningPolynomial(A) is not defined, the primitive
    element has components a and a+1), Cl(O) = Z/5 x Z/5, [O:E] ~ 2.6 10^16 = 2^16 3^8 11^2 23^2 941. E has 8 singular
    primes: three above 3, two above 23, one each above 2 (residue field F_4), 11 (residue field F_121) and 941 (split),
    so the lattice of overorders is the product of 8 local lattices. The biggest OverOrders computation of the suite:
    the Cartesian product, the orders generated from pairs of Z-bases and the final sorting dominate.
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
K:=NumberField(x^4+x^2+529);
A:=EtaleAlgebra([K,K]);
E:=EquationOrder(A);
t0:=Cputime();
    oo:=OverOrders(E);
t:=Cputime(t0);
assert #oo eq 36608;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
