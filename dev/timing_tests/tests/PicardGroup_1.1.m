/*
    PicardGroup_1.1: PicardGroup (GRH) of the equation order E of f = x^4-1000x^3-1000x^2-1000x-1000, computed on 5
    fresh copies of the algebra (the result is cached in the order), compared with Magma's PicardGroup of the RngOrd E.
    Why: f is irreducible of signature (2,1); the maximal order O has class group Z/3548, so the discrete logarithms
    in Pic(O) are not trivial, and #Pic(E) = 3548000. [O:E] = 1000; E has 2 singular primes, above 2 and 5, both
    totally ramified in O (e = 4, residue field F_p); E is Gorenstein. This was the slowest example of the package for
    a long time, with erratic timings caused by the random choices in CoprimeRepresentative; the RngOrd version is
    the natural yardstick.
*/
AttachSpec(spec);
SetSeed(1);
SetClassGroupBounds("GRH");
_<x>:=PolynomialRing(Integers());
f:=x^4-1000*x^3-1000*x^2-1000*x-1000;
K2:=NumberField(f);
E2:=EquationOrder(K2);
P2:=PicardGroup(E2);
assert #P2 eq 3548000;
t:=0;
for i in [1..5] do
    A:=EtaleAlgebra(f);
    E:=EquationOrder(A);
    t0:=Cputime();
        P,p:=PicardGroup(E : GRH:=true);
    t+:=Cputime(t0);
    assert #P eq #P2;
end for;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
