/*
    ICM_5.7: ICM (GRH) of the equation orders E = Z[F] of the 148 irreducible Weil polynomials x^4+ax^2+41^2 with a
    coprime to 41 (ordinary, q = 41).
    Why: a family of CM quartic fields with class groups of all sizes, from trivial to more than 100 elements, and
    equation orders of different indices; each ICM is small (from a few classes to 5628, 123279 in total), so the test
    measures the whole chain PicardGroup, WKICM, ICM on many different orders. Replaces ICM_5.2 (only
    x^4-17x^2+841, q = 29, 920 classes, about 1.5 s with 9362a5e).
*/
AttachSpec(spec);
SetSeed(1);
SetClassGroupBounds("GRH");
_<x>:=PolynomialRing(Integers());
q:=41;
polys:=[ x^4+a*x^2+q^2 : a in [-2*q+1..2*q-1] | Gcd(a,q) eq 1 and IsIrreducible(x^4+a*x^2+q^2) ];
nums:=[];
t:=0;
for f in polys do
    A:=EtaleAlgebra(f);
    E:=EquationOrder(A);
    t0:=Cputime();
        icm:=ICM(E : GRH:=true);
    t+:=Cputime(t0);
    Append(~nums,#icm);
end for;
assert #polys eq 148 and &+nums eq 123279 and Max(nums) eq 5628;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
