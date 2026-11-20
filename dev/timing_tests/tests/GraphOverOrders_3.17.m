/*
    GraphOverOrders_3.17: GraphOverOrders of the equation orders of four algebras of dimension 8 and 10, when the
    overorders are known but not their minimal inclusions, as after LoadWKICM: 42558 vertices and 216503 edges in
    total.
    Why: the Weil polynomials x^10-x^9+4x^8-6x^7+8x^6-16x^5+16x^4-24x^3+32x^2-16x+32 (q = 2, four components, 2
    singular primes, 2286 vertices) and x^8-2x^7+7x^6-14x^5+40x^4-56x^3+112x^2-128x+256 (q = 2, the order of
    WKICM_4.7, 3 singular primes, 16320 vertices), whose local lattices are large and deep, and the products of five
    and four imaginary quadratic fields (x^2+1)(x^2+2)(x^2+3)(x^2+5)(x^2+7) and (x^2+5)(x^2+7)(x^2+11)(x^2+13) (5 and
    4 singular primes, 14400 and 9552 vertices): the minimal inclusions are computed with MinimalOverOrdersAtPrime on
    every P-overorder. The overorders are computed in another copy of the algebra (not timed). Replaces
    GraphOverOrders_3.9 (the first order, after OverOrders, about 3 s with 9362a5e).
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
fs:=[ x^10-x^9+4*x^8-6*x^7+8*x^6-16*x^5+16*x^4-24*x^3+32*x^2-16*x+32,
      x^8-2*x^7+7*x^6-14*x^5+40*x^4-56*x^3+112*x^2-128*x+256, (x^2+1)*(x^2+2)*(x^2+3)*(x^2+5)*(x^2+7),
      (x^2+5)*(x^2+7)*(x^2+11)*(x^2+13) ];
// < the algebra and a fresh copy of it, number of overorders = vertices, number of edges >
nvs:=[ 2286, 16320, 14400, 9552 ];
nes:=[ 8985, 83966, 76880, 46672 ];
data:=[* <EtaleAlgebra(fs[i]),EtaleAlgebra(fs[i]),nvs[i],nes[i]> : i in [1..4] *];
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
