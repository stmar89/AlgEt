/*
    CohenMacaulay_9.4: CohenMacaulayType, IsGorenstein and IsBass of all the 9552 overorders of the equation order E
    of f = (x^2+5)(x^2+7)(x^2+11)(x^2+13) (OverOrders not timed).
    Why: four imaginary quadratic fields; [O:E] = 2^20 3^4 and E has 4 singular primes. The overorders have
    Cohen-Macaulay type up to 5 and 7768 of them are not Gorenstein: each type needs the trace dual ideal and the
    index of P*S^t in S^t at the non-Gorenstein primes, a test of TraceDualIdeal, ColonIdeal and Index over many
    small orders. Replaces CohenMacaulay_9.1 (the first three fields, 276 overorders, about 6 s with 9362a5e).
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=(x^2+5)*(x^2+7)*(x^2+11)*(x^2+13);
A:=EtaleAlgebra(f);
E:=EquationOrder(A);
oo:=OverOrders(E);
assert #oo eq 9552;
t0:=Cputime();
    cm:=[ CohenMacaulayType(T) : T in oo ];
    gor:=[ IsGorenstein(T) : T in oo ];
    bass:=[ IsBass(T) : T in oo ];
t:=Cputime(t0);
assert {* c : c in cm *} eq {* 1^^1784, 2^^4360, 3^^2768, 4^^592, 5^^48 *};
assert #[ g : g in gor | not g ] eq 7768 and #[ b : b in bass | b ] eq 320;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
