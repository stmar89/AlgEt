/*
    ICM_5.6: ICM (GRH) of the equation order E of f = (x^2+23)(x^2+92): 171432 ideal classes.
    Why: both components are Q(sqrt(-23)), whose class group is Z/3, given by different defining polynomials, so the
    package sees two distinct fields that it cannot identify; [O:E] = 2^3 3^2 23^2 and E has 5 singular primes.
    WKICM(E) has 72 classes and Pic(E) has 18216 elements, so the time goes into the products of the weak classes
    with the Picard representatives and the sorting of the classes by multiplicator ring. Replaces ICM_5.1
    ((x^2+5)(x^2+20), 1356 classes, about 3 s with 9362a5e).
*/
AttachSpec(spec);
SetSeed(1);
SetClassGroupBounds("GRH");
_<x>:=PolynomialRing(Integers());
f:=(x^2+23)*(x^2+92);
A:=EtaleAlgebra(f);
E:=EquationOrder(A);
t0:=Cputime();
    icm:=ICM(E : GRH:=true);
t:=Cputime(t0);
assert #icm eq 171432;
assert #WKICM(E) eq 72 and #PicardGroup(E : GRH:=true) eq 18216;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
