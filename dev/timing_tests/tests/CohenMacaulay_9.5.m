/*
    CohenMacaulay_9.5: CohenMacaulayType of all the 16320 overorders of the equation order E of
    f = x^8-2x^7+7x^6-14x^5+40x^4-56x^3+112x^2-128x+256, a Weil polynomial for q = 2 (the order of WKICM_4.7;
    OverOrders not timed).
    Why: [O:E] = 2^18 3^6 and E has 3 singular primes; 14250 of the overorders are not Gorenstein and the types go up
    to 5. Many more orders than in CohenMacaulay_9.4, in dimension 8: the per-order overhead of TraceDualIdeal,
    ColonIdeal and Index. Replaces CohenMacaulay_9.2 (the 2288 overorders of (x^4+16)(x^4+81), about 18 s with
    9362a5e).
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=x^8-2*x^7+7*x^6-14*x^5+40*x^4-56*x^3+112*x^2-128*x+256;
A:=EtaleAlgebra(f);
E:=EquationOrder(A);
oo:=OverOrders(E);
assert #oo eq 16320;
t0:=Cputime();
    cm:=[ CohenMacaulayType(T) : T in oo ];
t:=Cputime(t0);
assert Max(cm) eq 5 and #[ c : c in cm | c gt 1 ] eq 14250;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
