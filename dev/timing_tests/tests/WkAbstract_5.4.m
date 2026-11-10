/*
    WkAbstract_5.4: MultiplicationTable of the abstract weak equivalence class monoid W of R = Z[F,16/F] for
    f = x^6-3x^5-3x^4+65x^3-48x^2-768x+4096 (W itself is computed before, not timed): 204 classes, 20910 products.
    Why: A is the product of an imaginary quadratic and a quartic CM field (see OverOrders_3.2), [O:R] = 2^2 3^9 and R
    has 3 singular primes and 144 overorders, 60 of them non-Gorenstein. Every product in W multiplies the two
    representative ideals and then looks the result up among the classes with the same multiplicator ring with
    IsWeaklyEquivalent (two colon ideals and a product each): with 204 classes the table takes minutes (~270 s in
    2025), which makes it the benchmark for the cost of '*' and '!' in AlgEtQWECM.
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
f:=x^6 - 3*x^5 - 3*x^4 + 65*x^3 - 48*x^2 - 768*x + 4096;
A:=EtaleAlgebra(f);
F:=PrimitiveElement(A);
R:=Order([F,16/F]);
W,w:=WeakEquivalenceClassMonoidAbstract(R);
assert #W eq 204;
t0:=Cputime();
    T:=MultiplicationTable(W);
t:=Cputime(t0);
assert #Keys(T) eq 204*205 div 2;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
