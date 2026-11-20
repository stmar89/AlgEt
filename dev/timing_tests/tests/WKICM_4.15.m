/*
    WKICM_4.15: WKICM of R = Z[F,q/F] for 12 Weil polynomials of degree 8 and 10 with q = 3, 4, 5, 7 and 9, each
    after OverOrders(R : populateoo_in_oo) (not timed): 34529 weak equivalence classes in total.
    Why: R = Z[F,q/F] is the order that occurs for abelian varieties (it contains the Verschiebung q/F); these 12
    have between 2 and 9 singular primes, indices [O:R] from 2^4 3^6 11^2 to 2^22 3^4, Cohen-Macaulay types > 1 and
    from 147 to 6144 overorders, so the test runs the recursion over the overorders and the gluing of local classes
    on many different orders. Found by a search over the Weil polynomials of degree <= 10 and q <= 9 for orders
    Z[F,q/F] of big index. Replaces WKICM_4.1 (Z[F,2/F] for x^8+16, 31 classes, on 10 fresh copies, about 4 s with
    9362a5e).
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
// < f, q, number of weak equivalence classes of Z[F,q/F] >
data:=[
    <x^10-x^9+3*x^8-4*x^7+6*x^6-6*x^5+18*x^4-36*x^3+81*x^2-81*x+243,3,3648>,
    <x^10+x^8+9*x^6+27*x^4+27*x^2+243,3,3696>,
    <x^8+2*x^6+3*x^4+18*x^2+81,3,285>,
    <x^8+x^7+3*x^6+27*x^2+27*x+81,3,420>,
    <x^8+x^7+4*x^6+4*x^5+24*x^4+16*x^3+64*x^2+64*x+256,4,816>,
    <x^8+3*x^7+7*x^6+14*x^5+24*x^4+56*x^3+112*x^2+192*x+256,4,360>,
    <x^8-x^7-3*x^5-2*x^4-15*x^3-125*x+625,5,1432>,
    <x^10-x^8-4*x^7+13*x^6-4*x^5+65*x^4-100*x^3-125*x^2+3125,5,2816>,
    <x^8+2*x^6-45*x^4+98*x^2+2401,7,6840>,
    <x^8+2*x^7+3*x^6+16*x^5+64*x^4+112*x^3+147*x^2+686*x+2401,7,3200>,
    <x^8+2*x^7-12*x^5-63*x^4-108*x^3+1458*x+6561,9,4872>,
    <x^10+3*x^9+7*x^8+27*x^7+91*x^6+348*x^5+819*x^4+2187*x^3+5103*x^2+19683*x+59049,9,6144>
];
assert &+[ d[3] : d in data ] eq 34529;
t:=0;
for d in data do
    f,q,n:=Explode(d);
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    R:=Order([F,q/F]);
    oo:=OverOrders(R : populateoo_in_oo:=true);
    t0:=Cputime();
        wk:=WKICM(R);
    t+:=Cputime(t0);
    assert #wk eq n;
end for;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
