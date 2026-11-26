/*
    ZBasisLLL_2.5: ZBasisLLL of 20 ideals with a small index in R and huge Z-bases, 4 in each of five algebras of
    dimension 12 to 24, 100 rounds on fresh copies.
    Why: the stored basis is a*HNF(a^-1 L), where L is an ideal of small index of the equation order R and a = b^k
    is a power of an element b with small coordinates: exactly the shape of the bases built by SmallRepresentative
    (a times the HNF of an ideal), with entries of 800 to 6000 digits while the lattice L has a small index.
    ZBasisLLL has to notice that the HNF is the better starting point for the LLL. The algebras:
    (x^8+16)(x^6+8x^5+50x^4+200x^3+1250x^2+5000x+15625) (a CM octic and three imaginary quadratic fields, k = 640),
    (x^6+3)(x^6-2)(x^4+5) (k = 1280), the 72nd cyclotomic field (dimension 24, k = 640), x^12-x-1 (k = 1280) and
    (x^4+x^2+529)(x^5-x-1)(x^3-2)(x^2+5) (k = 640).
    Inputs (fixed, see README): L and b come from the Lehmer generator below (state 1 in each algebra); the bases are
    a times the HNF basis of a^-1 L, the same in every version. Replaces ZBasisLLL_2.2 (20 ideals in the first
    algebra, with entries of about 100 digits, 250 rounds, about 6 s with 9362a5e).
*/
AttachSpec(spec);
// Stores zb as the Z-basis of J without any computation: before 2026-11-26 the attribute ZBasis is the
// sequence, since then a record (zbasis_fmt in Ord.m), whose format is read from an existing record.
set_zbasis:=procedure(J,zb)
    if "CanonicalZBasis" in GetAttributes(AlgEtQIdl) then
        S:=Order(J); _:=ZBasis(OneIdeal(S));
        J`ZBasis:=rec< Format(OneIdeal(S)`ZBasis) | A:=Algebra(J), seq:=zb, lll:=false, hnf:=false >;
    else
        J`ZBasis:=zb;
    end if;
end procedure;
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
data:=[ <(x^8+16)*(x^6+8*x^5+50*x^4+200*x^3+1250*x^2+5000*x+15625),640>, <(x^6+3)*(x^6-2)*(x^4+5),1280>,
        <x^24-x^12+1,640>, <x^12-x-1,1280>, <(x^4+x^2+529)*(x^5-x-1)*(x^3-2)*(x^2+5),640> ];
inputs:=[* *];
chks:=[];
for d in data do
    f,k:=Explode(d);
    A:=EtaleAlgebra(f);
    R:=EquationOrder(A);
    s:=1;
    for i in [1..4] do
        repeat bs,s:=lehmer_elts(A,1,s); b:=bs[1]; until not IsZeroDivisor(b);
        gs,s:=lehmer_elts(A,2,s);
        L:=Ideal(R,gs);
        a:=b^k;
        zb:=[ a*z : z in ZBasis(Ideal(R,ZBasis((1/a)*L))) ]; // a times the HNF basis, the same in every version
        Append(~chks,checksum(zb));
        Append(~inputs,<R,zb,L>);
    end for;
end for;
printf "### INPUT %o\n",chks;
assert chks eq [ 764809953053649185, 566465500215078976, 431929235987359778, 418591247334423944,
                 686286758964178068, 193059848975645232, 1028854964610433754, 1376544078807916767,
                 706653502304138474, 2031581952532350986, 1103585459289029336, 1636988266156792491,
                 1369662899173263818, 1654465627066469576, 1564727684664170948, 2087567474725106677,
                 2027644092864496008, 461757376153429352, 1957409402532488378, 137724587596916147 ];
t:=0;
for r in [1..100] do
    for inp in inputs do
        R,zb,L:=Explode(inp);
        J:=Ideal(R,zb); // a fresh copy, since ZBasisLLL does nothing on a reduced basis
        set_zbasis(J,zb);
        t0:=Cputime();
            ZBasisLLL(J);
        t+:=Cputime(t0);
        if r eq 1 then assert Ideal(R,ZBasis(J)) eq L; end if; // checked once, to save time
    end for;
end for;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
