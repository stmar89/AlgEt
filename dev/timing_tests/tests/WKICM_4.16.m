/*
    WKICM_4.16: WKICM_bar of all the overorders of Cohen-Macaulay type >= 3 of the equation orders of three algebras:
    1118 orders, 11624 weak equivalence classes with the given multiplicator ring (OverOrders not timed).
    Why: an order of type >= 3 goes through the general branch of WKICM_bar: local computations at each singular
    prime (orders S+P^kO), then the gluing of the local classes. The algebras: the Weil polynomial
    x^10-x^9+4x^8-6x^7+8x^6-16x^5+16x^4-24x^3+32x^2-16x+32 (q = 2, four components; 702 orders of type 3 and 4),
    (x^4+16)(x^4+81) (two CM quartic fields; 368 orders of type 3) and (x^2+5)(x^2+7)(x^2+11) (three imaginary
    quadratic fields; 48 orders of type 3 and 4). The orders are the overorders of type >= 3, a set that does not
    depend on the code. Replaces WKICM_4.2 (WKICM_bar of one order of type 3, on 6 fresh copies, about 5 s with
    9362a5e).
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
// < f, number of overorders of type >= 3 of the equation order, total number of classes >
data:=[ <x^10-x^9+4*x^8-6*x^7+8*x^6-16*x^5+16*x^4-24*x^3+32*x^2-16*x+32,702,9016>, <(x^4+16)*(x^4+81),368,2192>,
        <(x^2+5)*(x^2+7)*(x^2+11),48,416> ];
t:=0;
for d in data do
    f,n,m:=Explode(d);
    A:=EtaleAlgebra(f);
    E:=EquationOrder(A);
    hi:=[ S : S in OverOrders(E) | CohenMacaulayType(S) ge 3 ];
    assert #hi eq n;
    t0:=Cputime();
        wk:=[ WKICM_bar(S) : S in hi ];
    t+:=Cputime(t0);
    assert &+[ #w : w in wk ] eq m;
end for;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
