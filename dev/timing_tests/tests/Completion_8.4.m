/*
    Completion_8.4: Uniformizers of all the primes of the maximal order above p <= 20 in the 72nd cyclotomic field
    (43 primes) and above p <= 30 in the algebra of (x^6+3)(x^6-2)(x^4+5) (69 primes), then Completion at each of
    them and the check of the valuations.
    Why: Uniformizers builds, for each prime, an element that is a uniformizer there and a unit at all the other
    primes (quotients P/P^2 and a CRT with all the moduli), so its cost grows with the number of primes. Q(zeta_72)
    has dimension 24, with 2 and 3 ramified (e = 4 and 6), residue degree 6 at 5, 7, 11, 13 and 2 at 17, 19; the second
    algebra has three components of dimension 6, 6 and 4, with primes of many ramification indices and residue
    degrees.
    Replaces Completion_8.2 (18 primes in the algebra of (x^8+16)(x^8+81), about 3 s with 9362a5e).
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
data:=[ <x^24-x^12+1,20,43>, <(x^6+3)*(x^6-2)*(x^4+5),30,69> ];
t:=0;
for d in data do
    f,b,n:=Explode(d);
    A:=EtaleAlgebra(f);
    O:=MaximalOrder(A);
    pp:=&cat[ PrimesAbove(p*O) : p in PrimesUpTo(b) ];
    assert #pp eq n;
    t0:=Cputime();
        unifs:=Uniformizers(pp);
        for iP->P in pp do
            AP,mP:=Completion(P);
            for it->u in unifs do
                assert Valuation(mP(u)) eq (iP eq it select 1 else 0);
            end for;
        end for;
    t+:=Cputime(t0);
end for;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
