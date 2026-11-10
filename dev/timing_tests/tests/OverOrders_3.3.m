/*
    OverOrders_3.3: OverOrders of the equation order E of f = x^6+8x^5+50x^4+200x^3+1250x^2+5000x+15625: 3200 overorders.
    Why: f is a Weil polynomial for q = 25 and factors as (x^2-5x+25)(x^2+5x+25)(x^2+8x+25): the algebra is the
    product of three imaginary quadratic fields, two of them isomorphic to Q(sqrt(-3)) (with different defining
    polynomials) and one to Q(i); Cl(O) is trivial. [O:E] ~ 1.8 10^11; E has 5 singular primes (above 2, 3, 5 and
    two above 13), so the lattice of overorders is a product of 5 small local lattices (compare with 3.2, three deep
    local lattices of about the same total size): the Cartesian product and the final sorting weigh more here.
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=x^6 + 8*x^5 + 50*x^4 + 200*x^3 + 1250*x^2 + 5000*x + 15625;
A:=EtaleAlgebra(f);
E:=EquationOrder(A);
t0:=Cputime();
    oo:=OverOrders(E);
t:=Cputime(t0);
assert #oo eq 3200;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
