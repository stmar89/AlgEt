/*
    WKICM_4.12: WKICM of the equation order E of f = x^4 + 3^5 5^3 7^2: 1500 classes, 624 overorders.
    Why: a CM quartic field with Cl(O) = Z/2 x Z/4; [O:E] = 2^3 3^6 5^3 7^2 and E has 4 singular primes with four
    different local pictures: above 2 (one ramified and one inert prime of O), above 3 and above 5 (totally ramified),
    above 7 (e = f = 2). Pic(E) has 17280 elements. A mid-size example where the local lattices are of different
    nature; its 4 singular primes also make IsPrincipal slow (see IsPrincipal_6.2).
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=x^4+3^5*5^3*7^2;
A:=EtaleAlgebra(f);
E:=EquationOrder(A);
t0:=Cputime();
    wk:=WKICM(E);
t:=Cputime(t0);
assert #wk eq 1500;
assert #OverOrders(E) eq 624;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
