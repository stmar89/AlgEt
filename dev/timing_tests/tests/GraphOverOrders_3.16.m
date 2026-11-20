/*
    GraphOverOrders_3.16: GraphOverOrders of the equation order of f = (x^4+16)(x^4+81)(x^2+5): 54912 vertices,
    374720 edges, when the overorders are known but not their minimal inclusions, as after LoadWKICM.
    Why: Q(zeta_8), Q(zeta_12) and Q(sqrt(-5)); [O:E] = 2^8 3^6 5^4 13^4 41^2 53^2 and E has 10 singular primes,
    so the graph is the product of 10 small local graphs: computing them is cheap, and the time goes into the
    matching of the 10-tuples of local overorders with the positions in OverOrders and into the 374720 edges. The
    overorders are computed in another copy of the algebra (not timed). Replaces GraphOverOrders_3.6 (5 singular
    primes, 3200 vertices, after OverOrders, about 1 s with 9362a5e).
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
f:=(x^4+16)*(x^4+81)*(x^2+5);
// < the algebra and a fresh copy of it, number of overorders = vertices, number of edges >
data:=[* <EtaleAlgebra(f),EtaleAlgebra(f),54912,374720> *];
t:=0;
for d in data do
    A1,A2,nvert,nedge:=Explode(d);
    oo:=OverOrders(EquationOrder(A1));
    R2:=copy_with_overorders(oo,A2);
    assert R2 eq EquationOrder(A2);
    t0:=Cputime();
        G:=GraphOverOrders(R2);
    t+:=Cputime(t0);
    assert #Vertices(G) eq nvert and #Edges(G) eq nedge;
end for;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
