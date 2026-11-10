/*
    GraphOverOrders_3.4: GraphOverOrders of the equation order of f = x^4-10000x^3-10000x^2-10000x-10000, after
    OverOrders (not timed): 297 vertices, 950 edges; computed on 12 fresh copies of the algebra.
    Why: 2 singular primes, so the graph is the product of two local graphs and the code has to match the Cartesian
    product of the local lattices with the sorted output of OverOrders. See OverOrders_3.1.
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=x^4-10000*x^3-10000*x^2-10000*x-10000;
t:=0;
for i in [1..12] do
    A:=EtaleAlgebra(f);
    E:=EquationOrder(A);
    oo:=OverOrders(E);
    t0:=Cputime();
        G:=GraphOverOrders(E);
    t+:=Cputime(t0);
    assert #Vertices(G) eq #oo and #Vertices(G) eq 297;
    assert #Edges(G) eq 950;
end for;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
