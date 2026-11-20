/*
    OverOrders_3.2: OverOrders of the equation order E of f = x^6-3x^5-3x^4+65x^3-48x^2-768x+4096: 3312 overorders.
    Why: f is a Weil polynomial for q = 16 and factors as (x^2-7x+16)(x^4+4x^3+9x^2+64x+256), so the algebra is the
    product of an imaginary quadratic and a quartic CM field. [O:E] = 322486272 = 2^19 3^5; E has 3 singular primes,
    two above 2 (residue field F_2; above each of them one unramified and one ramified prime of O) and one above 3
    (three ramified primes of O above it). Deep local lattices: the biggest lattice of the suite together with 3.3,
    and a good profile of MinimalOverOrdersAtPrime.
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=x^6 - 3*x^5 - 3*x^4 + 65*x^3 - 48*x^2 - 768*x + 4096;
A:=EtaleAlgebra(f);
E:=EquationOrder(A);
t0:=Cputime();
    oo:=OverOrders(E);
t:=Cputime(t0);
assert #oo eq 3312;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
