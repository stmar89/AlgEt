/*
    Completion_8.2: Uniformizers of the 18 primes of the maximal order of f = (x^8+16)(x^8+81) above 2, 3, 5 and 7,
    then Completion at each of them and the check of the valuations.
    Why: Uniformizers builds, for each prime, an element that is a uniformizer there and a unit at the other 17 primes
    (quotients P/P^2 and a CRT with 18 moduli), in an algebra of dimension 16 with 2 components; the primes above 2
    and 3 are ramified (e = 8 and e = 2), the others have residue degree 4.
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=(x^8+16)*(x^8+81);
A:=EtaleAlgebra(f);
O:=MaximalOrder(A);
pp:=PrimesAbove(2*3*5*7*O);
assert #pp eq 18;
t0:=Cputime();
    unifs:=Uniformizers(pp);
    for iP->P in pp do
        AP,mP:=Completion(P);
        for it->u in unifs do
            assert Valuation(mP(u)) eq (iP eq it select 1 else 0);
        end for;
    end for;
t:=Cputime(t0);
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
