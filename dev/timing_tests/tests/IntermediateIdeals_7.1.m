/*
    IntermediateIdeals_7.1: the fractional E-ideals between the conductor f and E O, for the equation order E of
    f = x^4+291x^3-988x^2-1000x-1000: all of them (2736), those with multiplicator ring E (120) and those with trivial
    extension to O (480).
    Why: a quartic of signature (2,1) with trivial Cl(O), [O:E] = 120 = 2^3 3 5 and 3 singular primes above 2, 3 and 5
    (2 has an inert and a split prime above it, 3 inert, 5 ramified), so E O / f is a module with several simple
    constituents. The three variants exercise the minimal-submodule search (MinimalSubmodules over F_p), the
    multiplicator-ring filter and the maximal-submodule search; the first is the slowest block of the fast tests.
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=x^4+291*x^3-988*x^2-1000*x-1000;
K:=EtaleAlgebra(f);
E:=EquationOrder(K);
O:=MaximalOrder(K);
EO:=E!!OneIdeal(O);
ff:=Conductor(E);
t0:=Cputime();
    nn:=IntermediateIdeals(EO,ff);
    n1:=IntermediateIdeals(EO,ff : PrescribedMultiplicatorRing:=true);
    n2:=IntermediateIdeals(EO,ff,O);
t:=Cputime(t0);
assert #nn eq 2736 and #n1 eq 120 and #n2 eq 480;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
