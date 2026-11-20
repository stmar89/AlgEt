/*
    OverOrders_3.12: OverOrders of the equation order E of f = (x^4+x^2+529)(x^4+16)(x^2+5): 180224 overorders.
    Why: three components, the CM quartic field of the Weil polynomial x^4+x^2+529 (q = 23, see OverOrders_3.7),
    Q(zeta_8) given by x^4+16 and Q(sqrt(-5)); [O:E] = 2^8 3^6 5^2 13^2 23 41^2 61^2 4049^2 and E has 13 singular
    primes, several of them coming from the resultants of the factors. The lattice of overorders is the product of 13
    small local lattices, so most of the time goes into the sums of local overorders and the final sorting.
    Replaces OverOrders_3.3 (5 singular primes, 3200 overorders, about 3 s with 9362a5e).
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=(x^4+x^2+529)*(x^4+16)*(x^2+5);
A:=EtaleAlgebra(f);
E:=EquationOrder(A);
t0:=Cputime();
    oo:=OverOrders(E);
t:=Cputime(t0);
assert #oo eq 180224;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
