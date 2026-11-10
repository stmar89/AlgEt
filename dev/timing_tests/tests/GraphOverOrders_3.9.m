/*
    GraphOverOrders_3.9: GraphOverOrders of the equation order E of f = x^10-x^9+4x^8-6x^7+8x^6-16x^5+16x^4-24x^3+32x^2
    -16x+32, a Weil polynomial for q = 2 of degree 10, after OverOrders (not timed): 2286 vertices, 8985 edges.
    Why: f = (x^2-x+2)(x^2+2)(x^2+2x+2)(x^4-2x^3+2x^2-4x+4), so A has 4 components (Q(sqrt(-7)), Q(sqrt(-2)), Q(i) and
    a quartic CM field); Cl(O) is trivial, [O:E] = 2^17 3^6, and E has 2 singular primes: one above 3 with residue
    field F_9 and one above 2 with residue field F_2 (below one unramified and three ramified primes of O). The lattice
    of overorders is the product of two local lattices, one of them large; building the graph from the local graphs
    takes 4 times longer than OverOrders itself (88 s vs 21 s in 2025).
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=x^10 - x^9 + 4*x^8 - 6*x^7 + 8*x^6 - 16*x^5 + 16*x^4 - 24*x^3 + 32*x^2 - 16*x + 32;
A:=EtaleAlgebra(f);
E:=EquationOrder(A);
oo:=OverOrders(E);
assert #oo eq 2286;
t0:=Cputime();
    G:=GraphOverOrders(E);
t:=Cputime(t0);
assert #Vertices(G) eq 2286 and #Edges(G) eq 8985;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
