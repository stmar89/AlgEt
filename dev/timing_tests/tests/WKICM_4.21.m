/*
    WKICM_4.21: WKICM of the equation orders E of f = x^4 + 3^5 5^6 7^4, x^4 + 3^9 5^6 7^2 and
    x^4 + 2^3 3^5 5^6 7^2 11^2: 102000, 79650 and 57600 classes.
    Why: CM quartic fields with small class groups (Z/2, (Z/2)^3 and (Z/2)^5 x Z/4); [O:E] = 2^2 3^6 5^8 7^6,
    2^2 3^12 5^8 7^2 and 2^3 3^6 5^8 7^2 11^2, with 4, 4 and 5 singular primes of different local nature, and from
    19656 to 24570 overorders. Mid-size examples where the gluing of the local classes dominates. Replaces
    WKICM_4.12 (x^4 + 3^5 5^3 7^2, 1500 classes, about 2 s with 9362a5e).
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
data:=[ <x^4+3^5*5^6*7^4,102000>, <x^4+3^9*5^6*7^2,79650>, <x^4+2^3*3^5*5^6*7^2*11^2,57600> ];
t:=0;
for d in data do
    f,n:=Explode(d);
    A:=EtaleAlgebra(f);
    E:=EquationOrder(A);
    t0:=Cputime();
        wk:=WKICM(E);
    t+:=Cputime(t0);
    assert #wk eq n;
end for;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
