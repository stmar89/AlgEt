/*
    OverOrders_3.1: OverOrders of the equation order E of f = x^4-10000x^3-10000x^2-10000x-10000: 297 overorders,
    computed on 10 fresh copies of the algebra (the result is cached in the order).
    Why: a quartic of signature (2,1) with Cl(O) = Z/9764 and [O:E] = 10^6; E has 2 singular primes, above 2
    (totally ramified in O) and above 5 (which splits completely in O), so the lattice of overorders is the product
    of two local lattices of different nature. A quick test of the Hofmann-Sircana algorithm.
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=x^4-10000*x^3-10000*x^2-10000*x-10000;
t:=0;
for i in [1..10] do
    A:=EtaleAlgebra(f);
    E:=EquationOrder(A);
    t0:=Cputime();
        oo:=OverOrders(E);
    t+:=Cputime(t0);
    assert #oo eq 297;
end for;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
