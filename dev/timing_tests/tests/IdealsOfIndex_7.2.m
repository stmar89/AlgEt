/*
    IdealsOfIndex_7.2: IdealsOfIndex(E,N) for N = 1..30, for the equation order E of f = (x^8+16)(x^8+81), and the
    comparison of the fast and the slow method for the N coprime to the index.
    Why: [O:E] ~ 2.8 10^24 is divisible by 2, 3, 5 and 13, so for most N the index is not coprime to the conductor and
    the slow method (intermediate ideals of E/N E) is used; for the other N the ideals are built from the ideals of
    norm N of the two maximal orders. An algebra of dimension 16 with 2 components, where E/N E is large.
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=(x^8+16)*(x^8+81);
A:=EtaleAlgebra(f);
E:=EquationOrder(A);
O:=MaximalOrder(A);
ind:=Index(O,E);
t0:=Cputime();
    for N in [1..30] do
        ids:=IdealsOfIndex(E,N);
        if IsCoprime(N,ind) then
            assert Seqset(ids) eq Seqset(IdealsOfIndex(E,N : Method:="Slow"));
        end if;
    end for;
t:=Cputime(t0);
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
