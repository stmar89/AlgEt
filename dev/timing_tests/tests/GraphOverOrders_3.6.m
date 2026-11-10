/*
    GraphOverOrders_3.6: GraphOverOrders of the equation order of f = x^6+8x^5+50x^4+200x^3+1250x^2+5000x+15625,
    after OverOrders (not timed): 3200 vertices, 17376 edges.
    Why: 5 singular primes, so the 3200 overorders are a Cartesian product of 5 local lattices: the graph is built
    from the 5 local graphs and the matching of tuples with the sorted output of OverOrders is the expensive part.
    See OverOrders_3.3.
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=x^6 + 8*x^5 + 50*x^4 + 200*x^3 + 1250*x^2 + 5000*x + 15625;
A:=EtaleAlgebra(f);
E:=EquationOrder(A);
oo:=OverOrders(E);
t0:=Cputime();
    G:=GraphOverOrders(E);
t:=Cputime(t0);
assert #Vertices(G) eq #oo and #Vertices(G) eq 3200;
assert #Edges(G) eq 17376;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
