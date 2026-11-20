/*
    OverOrders_3.10: OverOrders of the equation orders of the 160 polynomials x^n - a k (x^(n-1)+...+x+1) for
    (n,a) = (3,10^6), (4,10^6), (5,10^4), (6,10^3) and k = 1..40 (all irreducible): 182494 overorders in total.
    Why: a family of fields of degree 3 to 6 and of different signatures; the singular primes of E are the primes
    dividing a k, at most 2 for most k, with different local pictures (at 2 and 5 the polynomial is far from
    Eisenstein, at the other primes dividing k it is often Eisenstein), so the lattices of overorders are products
    of two or three local lattices of different nature. A test of the Hofmann-Sircana algorithm on many small and
    medium orders. Replaces OverOrders_3.1 (x^4-10000(x^3+x^2+x+1), 297 overorders, on 10 fresh copies, about 3 s
    with 9362a5e).
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
nums:=[];
t:=0;
for na in [ <3,10^6>, <4,10^6>, <5,10^4>, <6,10^3> ] do
    n,a:=Explode(na);
    for k in [1..40] do
        f:=x^n-a*k*&+[ x^i : i in [0..n-1] ];
        if not IsIrreducible(f) then continue; end if;
        A:=EtaleAlgebra(f);
        E:=EquationOrder(A);
        t0:=Cputime();
            oo:=OverOrders(E);
        t+:=Cputime(t0);
        Append(~nums,#oo);
    end for;
end for;
assert #nums eq 160 and &+nums eq 182494;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
