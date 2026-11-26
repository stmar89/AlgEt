/*
    SmallRepresentative_2.4: SmallRepresentative of P^n, where P is a prime of the maximal order O of degree 1 whose
    class in Pic(O) has order n, in four fields of degree 2, 3, 4 and 6 with big class groups (one fresh copy of P^n
    each; the result is cached in the ideal).
    Why: P^n is principal, and it is exactly the kind of ideal that PicardGroup tests for principality to glue Pic(O)
    to the Picard group of a non-maximal order. Its Z-basis has entries with tens or hundreds of thousands of digits;
    SmallRepresentative multiplies by a short element of (O:P^n) and LLL-reduces the HNF of the result. The fields:
    Q(sqrt(-(10^10+19))) (Cl = Z/39809), x^3-10^7(x^2+x+1) (Cl = Z/465144), x^4-10^6(x^3+x^2+x+1) (Cl = Z/2 x Z/2 x
    Z/4 x Z/661864) and x^6-300(x^5+x^4+x^3+x^2+x+1) (Cl = Z/6 x Z/6 x Z/24090).
    Inputs (fixed, see README): P = pO + (F-c)O, with F the primitive element and c a root of f modulo p, p not
    dividing [O:Z[F]]; the input of SmallRepresentative is the HNF basis of P^n, the same in every version.
    Replaces SmallRepresentative_2.1 (the 3548th power of a prime in the quartic field of x^4-1000(x^3+x^2+x+1),
    1693 digits, about 1 s with 9362a5e).
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
// < f, p, c, the order n of the class of P = pO + (F-c)O >
data:=[ <x^2+10^10+19,5,4,39809>, <x^3-10^7*(x^2+x+1),23,4,465144>, <x^4-10^6*(x^3+x^2+x+1),3,1,661864>,
        <x^6-300*(x^5+x^4+x^3+x^2+x+1),17,7,24090> ];
chks:=[];
t:=0;
for d in data do
    f,p,c,n:=Explode(d);
    A:=EtaleAlgebra(f);
    O:=MaximalOrder(A);
    F:=PrimitiveElement(A);
    assert Index(O,EquationOrder(A)) mod p ne 0 and Evaluate(f,c) mod p eq 0;
    P:=Ideal(O,[ A!p, F-c ]);
    GO,gO:=PicardGroup(O : GRH:=true);
    assert Order(P@@gO) eq n;
    zb:=ZBasis(Ideal(O,ZBasis(P^n))); // the HNF basis, the same in every version
    Append(~chks,checksum(zb));
    J:=Ideal(O,zb); // a fresh copy
    set_zbasis(J,zb);
    t0:=Cputime();
        aJ,a:=SmallRepresentative(J);
    t+:=Cputime(t0);
    assert aJ subset OneIdeal(O) and aJ eq a*J;
end for;
printf "### INPUT %o\n",chks;
assert chks eq [ 144873859945431101, 446557048234135547, 882373066596285720, 469380492412655495 ];
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
