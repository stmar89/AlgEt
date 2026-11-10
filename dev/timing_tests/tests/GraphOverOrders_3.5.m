/*
    GraphOverOrders_3.5: GraphOverOrders of the equation order of f = x^6-3x^5-3x^4+65x^3-48x^2-768x+4096, after
    OverOrders (not timed): 3312 vertices, 12812 edges.
    Why: 3 singular primes with deep local lattices; with 3312 overorders the quadratic parts of the code (pairs of
    overorders, lookups in sequences) dominate. See OverOrders_3.2.
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=x^6 - 3*x^5 - 3*x^4 + 65*x^3 - 48*x^2 - 768*x + 4096;
A:=EtaleAlgebra(f);
E:=EquationOrder(A);
oo:=OverOrders(E);
t0:=Cputime();
    G:=GraphOverOrders(E);
t:=Cputime(t0);
assert #Vertices(G) eq #oo and #Vertices(G) eq 3312;
assert #Edges(G) eq 12812;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
