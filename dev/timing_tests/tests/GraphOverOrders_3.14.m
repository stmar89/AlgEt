/*
    GraphOverOrders_3.14: GraphOverOrders of the equation orders of the polynomials x^n - a k (x^(n-1)+...+x+1) for
    (n,a) = (3,10^6), (4,10^6), (5,10^4), (6,10^3) and k = 1..25 (96663 vertices and 378069 edges in total), when
    the overorders are known but not their minimal inclusions, as after LoadWKICM.
    Why: the family of OverOrders_3.10: 100 orders with two or three singular primes each, so each graph is the
    product of two or three local graphs; here GraphOverOrders has to compute the local graphs, with
    MinimalOverOrdersAtPrime on every P-overorder, and to match the tuples of local overorders with the positions in
    OverOrders. The overorders are computed in another copy of the algebra (not timed). Replaces
    GraphOverOrders_3.4 (one quartic of the family, after OverOrders, on 12 fresh copies, about 1 s with 9362a5e).
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
// Rebuilds the overorders oo of an order in a fresh copy A2 of its algebra, from their Z-bases, and assigns them to the
// copy R2 of the order, as LoadWKICM does: the overorders of R2 are known, but not their minimal inclusions, which
// GraphOverOrders(R2) has to compute (a fresh algebra is needed, since Order returns the known orders of A2).
copy_with_overorders:=function(oo,A2)
    zbs:=[ [ < Eltseq(c) : c in Components(z) > : z in ZBasis(S) ] : S in oo ];
    oo2:=[ Order([ A2!t : t in zb ] : Check:=0) : zb in zbs ];
    R2:=oo2[1];
    R2`OverOrders:=oo2;
    return R2;
end function;
nvert:=0; nedge:=0;
t:=0;
for na in [ <3,10^6>, <4,10^6>, <5,10^4>, <6,10^3> ] do
    n,a:=Explode(na);
    for k in [1..25] do
        f:=x^n-a*k*&+[ x^i : i in [0..n-1] ];
        if not IsIrreducible(f) then continue; end if;
        oo:=OverOrders(EquationOrder(EtaleAlgebra(f)));
        A2:=EtaleAlgebra(f);
        R2:=copy_with_overorders(oo,A2);
        assert R2 eq EquationOrder(A2);
        t0:=Cputime();
            G:=GraphOverOrders(R2);
        t+:=Cputime(t0);
        assert #Vertices(G) eq #oo;
        nvert+:=#oo; nedge+:=#Edges(G);
    end for;
end for;
assert nvert eq 96663 and nedge eq 378069;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
