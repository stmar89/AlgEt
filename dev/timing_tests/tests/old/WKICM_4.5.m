/*
    WKICM_4.5: WKICM of R = Z[F,25/F] for f = x^6-3x^5+30x^4-175x^3+750x^2-1875x+15625, a Weil polynomial for
    q = 25: 315 classes.
    Why: f = (x^2+5x+25)(x^4-8x^3+45x^2-200x+625), a product of Q(sqrt(-3)) and a quartic CM field with
    Cl(O) = Z/2 x Z/2. [O:R] = 360000 = 2^6 3^2 5^4; R has 3 singular primes (above 3, above 2 with residue field F_4,
    above 5), so the local-global gluing of WKICM runs over three local orders and the 315 classes are the product
    of the three sets of local classes.
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=x^6 - 3*x^5 + 30*x^4 - 175*x^3 + 750*x^2 - 1875*x + 15625;
A:=EtaleAlgebra(f);
F:=PrimitiveElement(A);
q:=Integers() ! Round(ConstantCoefficient(f)^(2/Degree(f)));
R:=Order([F,q/F]);
t0:=Cputime();
    wk:=WKICM(R);
t:=Cputime(t0);
assert #wk eq 315;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
