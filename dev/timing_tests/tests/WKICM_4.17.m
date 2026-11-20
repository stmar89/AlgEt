/*
    WKICM_4.17: WKICM of the equation orders E of x^9+3^5, x^6+2^7, x^8+2^5, x^9+3^4 and x^7+2^5: 2249 classes in
    total.
    Why: in each case the Newton polygon at p = 3 or 2 has one segment of slope coprime to the degree, so p is
    totally ramified in O and E has a unique singular prime, above p, with residue field F_p; [O:E] = 3^16, 2^15,
    2^14, 3^12, 2^12. WKICM uses the recursion over the minimal overorders, and many overorders have Cohen-Macaulay
    type > 2. The ramification is wild in all the cases (p divides the degree) except x^7+2^5. x^9+3^5 (1177 classes)
    takes most of the time. Replaces WKICM_4.3 (x^8+16, 173 classes, about 4 s with 9362a5e).
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
data:=[ <x^9+3^5,1177>, <x^6+2^7,491>, <x^8+2^5,263>, <x^9+3^4,177>, <x^7+2^5,141> ];
assert &+[ d[2] : d in data ] eq 2249;
t:=0;
for d in data do
    f,n:=Explode(d);
    A:=EtaleAlgebra(f);
    E:=EquationOrder(A);
    t0:=Cputime();
        wk:=WKICM(E);
    t+:=Cputime(t0);
    assert #wk eq n;
end for;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
