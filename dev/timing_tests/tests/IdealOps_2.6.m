/*
    IdealOps_2.6: arithmetic of ideals of the equation orders of five algebras of dimension 12 to 24: in each, 2500
    sums, 2500 products (with their Z-bases), 2500 colon ideals and 2500 intersections of pairs of ideals with 3
    generators, then the 2500 trace dual ideals.
    Why: the elementary operations behind every other test, on HNFs over Z of size 12 to 24 (48 x 24 for the sums of
    ideals with 3 generators). The algebras: x^12-x-1 (an irreducible polynomial with a big Galois group), the 72nd
    cyclotomic field (dimension 24), (x^8+16)(x^8+81) (two CM octic fields), (x^6+3)(x^6-2)(x^4+5) and
    (x^4+x^2+529)(x^5-x-1)(x^3-2)(x^2+5) (three and four components, signatures from totally complex to mixed).
    Inputs (fixed, see README): the generators come from the Lehmer generator below (state 1 in each algebra).
    Replaces IdealOps_2.3 (600 of each in (x^8+16)(x^8+81), products without Z-bases, about 6 s with 9362a5e).
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
// Fixed inputs (see README): a checksum of the coordinates of a sequence of elements, printed and asserted below.
checksum:=function(elts)
    s:=0;
    for z in elts do
        for c in Components(z) do
            for q in Eltseq(c) do
                s:=(1000003*s+7*Numerator(q)+Denominator(q)) mod (2^61-1);
            end for;
        end for;
    end for;
    return s;
end function;
// N pseudo-random elements with coefficients in [-3..3] on the powers of the primitive element of A, from a Lehmer
// generator with state s: the same elements in every Magma version, unlike Random. Returns them and the new state.
lehmer_elts:=function(A,N,s)
    F:=PrimitiveElement(A);
    pw:=[ F^i : i in [0..Dimension(A)-1] ];
    out:=[];
    for k in [1..N] do
        cs:=[];
        for i in [1..#pw] do
            s:=(48271*s) mod (2^31-1);
            Append(~cs,(s mod 7)-3);
        end for;
        Append(~out,&+[ cs[i]*pw[i] : i in [1..#pw] ]);
    end for;
    return out,s;
end function;
polys:=[ x^12-x-1, x^24-x^12+1, (x^8+16)*(x^8+81), (x^6+3)*(x^6-2)*(x^4+5), (x^4+x^2+529)*(x^5-x-1)*(x^3-2)*(x^2+5) ];
N:=2500;
chks:=[];
t:=0;
for f in polys do
    A:=EtaleAlgebra(f);
    E:=EquationOrder(A);
    s:=1;
    gens:=[];
    for i in [1..2*N] do
        gs,s:=lehmer_elts(A,3,s);
        Append(~gens,gs);
    end for;
    Append(~chks,checksum(&cat gens));
    ids:=[ Ideal(E,gens[i]) : i in [1..N] ];
    jds:=[ Ideal(E,gens[i]) : i in [N+1..2*N] ];
    t0:=Cputime();
        sums:=[ ids[i]+jds[i] : i in [1..N] ];
        prods:=[ ids[i]*jds[i] : i in [1..N] ];
        zbs:=[ ZBasis(I) : I in prods ];
        cols:=[ ColonIdeal(ids[i],jds[i]) : i in [1..N] ];
        meets:=[ ids[i] meet jds[i] : i in [1..N] ];
        duals:=[ TraceDualIdeal(I) : I in ids ];
    t+:=Cputime(t0);
    assert forall{ i : i in [1..20] | meets[i] subset ids[i] and ids[i] subset sums[i] and prods[i] subset meets[i] };
end for;
printf "### INPUT %o\n",chks;
assert chks eq [ 1080433119083292952, 208933541838104576, 178544711545183830, 1281122208179470344, 545696133064076970 ];
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
