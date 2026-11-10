/*
    ICM_5.2: ICM (GRH) of the equation order E of f = x^4-17x^2+841, a Weil polynomial for q = 29: 920 ideal classes.
    Why: a CM quartic field with Cl(O) = Z/8, the largest class group among the quartic Weil examples; [O:E] = 5^2 29
    and the 3 singular primes (two above 5, one above 29) all lie below a single inert prime of O (f = 2). E is Bass,
    so the 8 overorders give the 8 weak classes, and Pic(E) has 576 elements; the ideal class monoid has 920 classes
    (~10 s in 2025).
*/
AttachSpec(spec);
SetSeed(1);
SetClassGroupBounds("GRH");
_<x>:=PolynomialRing(Integers());
f:=x^4-17*x^2+841;
A:=EtaleAlgebra(f);
E:=EquationOrder(A);
t0:=Cputime();
    icm:=ICM(E : GRH:=true);
t:=Cputime(t0);
assert #icm eq 920;
assert #WKICM(E) eq 8 and #PicardGroup(E) eq 576;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
