/*
    IsPrincipal_6.3: IsPrincipal (GRH) on 11000 non-principal ideals of the equation orders of f = x^8 + 3^8 7 (5000),
    of A = K x K x K (4000) and of A = K x K (2000), K the CM quartic field of x^4+x^2+529 (PicardGroup not timed).
    Why: a non-principal ideal is recognized without computing a generator, from its residue in (O/f)^* modulo the
    images of E^* and O^*: 10 to 20 milliseconds each. The orders: x^8 + 3^8 7 (trivial Cl(O), unit rank 3,
    [O:E] = 2^6 3^28, #Pic(E) = 1259712), K x K and K x K x K (Cl(O) = (Z/5)^2 and (Z/5)^3, #Pic(E) = 382205952000
    and 67852818829292596730265600000, see PicardGroup_1.11). The ideals are products P^e Q of primes of degree 1,
    so their norms are coprime to the conductor and their classes are spread over Pic(E). IsPrincipal_6.4 has the
    principal case.
    Inputs (fixed, see README): the primes P = lE + (F-c)E for the primes l <= 400 not dividing [O:E] and the roots c
    of f modulo l, with F the primitive element; the indices of P and Q and the exponent e in [1..3] come from the
    Lehmer generator below (state 1); the products are taken through their HNF bases. Replaces IsPrincipal_6.1 (20
    non-principal ideals of the first order, about 1 s with 9362a5e).
*/
AttachSpec(spec);
SetSeed(1);
SetClassGroupBounds("GRH");
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
K:=NumberField(x^4+x^2+529);
data:=[* <[NumberField(x^8+3^8*7)],5000>, <[K,K,K],4000>, <[K,K],2000> *];
chks:=[];
t:=0;
for d in data do
    As,N:=Explode(d);
    A:=EtaleAlgebra(As);
    E:=EquationOrder(A);
    F:=PrimitiveElement(A);
    f:=PolynomialRing(Integers())!MinimalPolynomial(F);
    ind:=Index(MaximalOrder(A),E);
    prs:=[ Ideal(E,[A!l,F-c]) : l in PrimesUpTo(400), c in [0..399]
                                | c lt l and ind mod l ne 0 and Evaluate(f,c) mod l eq 0 ];
    s:=1;
    Js:=[];
    for i in [1..N] do
        s:=(48271*s) mod (2^31-1); a:=1+(s mod #prs);
        s:=(48271*s) mod (2^31-1); b:=1+(s mod #prs);
        s:=(48271*s) mod (2^31-1); e:=1+(s mod 3);
        Append(~Js,Ideal(E,ZBasis(Ideal(E,ZBasis(prs[a]^e*prs[b])))));
    end for;
    Append(~chks,checksum(&cat[ ZBasis(J) : J in Js ]));
    P,p:=PicardGroup(E : GRH:=true);
    t0:=Cputime();
        res:=[ IsPrincipal(J : GRH:=true) : J in Js ];
    t+:=Cputime(t0);
    assert forall{ r : r in res | not r };
end for;
printf "### INPUT %o\n",chks;
assert chks eq [ 1375992179323302312, 1028121379890966525, 1389009322186175501 ];
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
