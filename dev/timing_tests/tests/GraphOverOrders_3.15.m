/*
    GraphOverOrders_3.15: GraphOverOrders of the equation orders of five algebras, when the overorders are known but
    not their minimal inclusions, as after LoadWKICM: 89008 vertices and 501306 edges in total.
    Why: GraphOverOrders combines the graphs of the local lattices at the singular primes; here it has to compute
    them, with MinimalOverOrdersAtPrime on every P-overorder, and then match the tuples of local overorders with the
    positions in OverOrders. The algebras: A = K x K, K the CM quartic field of x^4+x^2+529 (8 singular primes, 36608
    vertices, see OverOrders_3.7), the Weil polynomial x^6-3x^5-3x^4+65x^3-48x^2-768x+4096 (q = 16, 3 singular primes
    with deep local lattices, 3312 vertices), (x^4+16)(x^4+81) (6 singular primes, 2288 vertices, see WKICM_4.8) and
    the CM quartic fields of x^4 + 3^5 5^6 7^4 and x^4 + 3^9 5^6 7^2 (4 singular primes, 22230 and 24570 vertices,
    see WKICM_4.21). The overorders are computed in another copy of the algebra (not timed). Replaces
    GraphOverOrders_3.5 (the sextic, after OverOrders, which leaves the minimal inclusions stored, about 1 s with
    9362a5e).
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
K:=NumberField(x^4+x^2+529);
f6:=x^6-3*x^5-3*x^4+65*x^3-48*x^2-768*x+4096;
// < the algebra and a fresh copy of it, number of overorders = vertices, number of edges >
f8:=(x^4+16)*(x^4+81);
g1:=x^4+3^5*5^6*7^4;
g2:=x^4+3^9*5^6*7^2;
data:=[* <EtaleAlgebra([K,K]),EtaleAlgebra([K,K]),36608,233728>, <EtaleAlgebra(f6),EtaleAlgebra(f6),3312,12812>,
        <EtaleAlgebra(f8),EtaleAlgebra(f8),2288,10656>, <EtaleAlgebra(g1),EtaleAlgebra(g1),22230,117750>,
        <EtaleAlgebra(g2),EtaleAlgebra(g2),24570,126360> *];
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
