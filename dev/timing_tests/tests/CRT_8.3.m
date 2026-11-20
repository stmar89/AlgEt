/*
    CRT_8.3: ChineseRemainderTheorem for 360 tuples of elements modulo the 20 primes of the conductor of the equation
    order E of f = (x^8+16)(x^8+81)(x^8+625) and their squares and fourth powers (120 tuples each), and for 50 pairs
    of elements modulo each of the 21 pairs of distinct conductor primes above the same rational prime.
    Why: the order of PicardGroup_1.9, of dimension 24, whose conductor has 20 primes above 2, 3, 5, 7, 13, 17 and 29,
    with residue fields from F_2 to F_29^4.
    The two-ideal CRT works with the HNF transformation of the stacked Z-bases (the integer shortcut does not apply
    to two primes above the same rational prime); the sequence version recurses with the intersections of the
    remaining ideals, which are computed again at every call: with 20 moduli and their powers this dominates.
    Inputs (fixed, see README): the primes are sorted by minimal integer, index and HNF basis, which does not depend
    on the code, and their HNF bases are checked; the elements come from the Lehmer generator below (state 1).
    Replaces CRT_8.1 (the 6 primes of the conductor of (x^8+16)(x^8+81), about 5 s with 9362a5e).
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
f:=(x^8+16)*(x^8+81)*(x^8+625);
A:=EtaleAlgebra(f);
assert [ DefiningPolynomial(K) : K in Components(A) ] eq [ x^8+16, x^8+81, x^8+625 ];
E:=EquationOrder(A);
pp:=PrimesAbove(Conductor(E));
// a key of each prime that does not depend on the code: minimal integer, index and the coordinates of the HNF basis
hnf:=func< P | ZBasis(Ideal(E,ZBasis(P))) >;
key:=func< P | [ MinimalInteger(P), Index(E,P) ] cat &cat[ &cat[ Eltseq(c) : c in Components(z) ] : z in hnf(P) ] >;
cmp:=function(P,Q)
    kP:=key(P); kQ:=key(Q);
    for i in [1..Min(#kP,#kQ)] do
        if kP[i] ne kQ[i] then return kP[i] lt kQ[i] select -1 else 1; end if;
    end for;
    return #kP-#kQ;
end function;
Sort(~pp,cmp);
assert [ MinimalInteger(P) : P in pp ] eq [ 2, 2, 3, 3, 3, 5, 5, 5, 7, 7, 7, 7, 13, 13, 17, 17, 17, 17, 29, 29 ];
assert [ Index(E,P) : P in pp ] eq [ 2, 2, 3, 81, 81, 5, 625, 625, 49, 49, 49, 49, 28561, 28561, 289, 289, 289, 289,
                                    707281, 707281 ];
chk_primes:=checksum(&cat[ hnf(P) : P in pp ]);
s:=1;
tuples:=[];
for i in [1..360] do
    es,s:=lehmer_elts(A,#pp,s);
    Append(~tuples,es);
end for;
pairs_idx:=[ [i,j] : i,j in [1..#pp] | i lt j and MinimalInteger(pp[i]) eq MinimalInteger(pp[j]) ];
pairs:=[];
for ij in pairs_idx do
    for k in [1..50] do
        repeat es,s:=lehmer_elts(A,1,s); a:=es[1]; until not a in pp[ij[1]];
        repeat es,s:=lehmer_elts(A,1,s); b:=es[1]; until not b in pp[ij[2]];
        Append(~pairs,<ij[1],ij[2],a,b>);
    end for;
end for;
chk:=checksum(&cat tuples cat &cat[ [ pr[3],pr[4] ] : pr in pairs ]);
printf "### INPUT %o\n",[ chk_primes, chk ];
assert [ chk_primes, chk ] eq [ 1459794737172388857, 1286116852685290967 ];
mods:=[ pp, [ P^2 : P in pp ], [ P^4 : P in pp ] ];
t0:=Cputime();
    for pr in pairs do
        e:=ChineseRemainderTheorem(pp[pr[1]],pp[pr[2]],pr[3],pr[4]);
    end for;
    for k in [1..3] do
        for tup in tuples[120*(k-1)+1..120*k] do
            e:=ChineseRemainderTheorem(mods[k],tup);
        end for;
    end for;
t:=Cputime(t0);
assert forall{ i : i in [1..#pp] | e-tuples[360][i] in mods[3][i] };
assert #pairs_idx eq 21;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
