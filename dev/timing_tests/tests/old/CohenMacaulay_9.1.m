/*
    CohenMacaulay_9.1: CohenMacaulayType, IsGorenstein and IsBass of all the 276 overorders of the equation order E of
    f = (x^2+5)(x^2+7)(x^2+11) (OverOrders not timed).
    Why: three imaginary quadratic fields, Cl(O) = Z/2, [O:E] = 2^10 3^2 and 3 singular primes (two above 3, which
    splits in Q(sqrt(-5)) and Q(sqrt(-7)), one above 2, below one ramified, two split and one inert prime of O); the
    overorders have Cohen-Macaulay type up to 4 and 172 of them are not Gorenstein. Each type needs the trace dual
    ideal and the index of P*S^t in S^t at the non-Gorenstein primes: a test of TraceDualIdeal, ColonIdeal and Index
    over many small orders.
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=(x^2+5)*(x^2+7)*(x^2+11);
A:=EtaleAlgebra(f);
E:=EquationOrder(A);
oo:=OverOrders(E);
assert #oo eq 276;
t0:=Cputime();
    cm:=[ CohenMacaulayType(T) : T in oo ];
    gor:=[ IsGorenstein(T) : T in oo ];
    bass:=[ IsBass(T) : T in oo ];
t:=Cputime(t0);
assert Max(cm) eq 4 and #[ g : g in gor | not g ] eq 172;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
