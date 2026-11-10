/*
    SmallRepresentative_2.1: SmallRepresentative of an ideal with a huge Z-basis, 5 times on fresh copies (the
    result is cached in the ideal).
    Why: the ideal is cI^n, where cI is a representative of a generator of Pic(O) coprime to the conductor and n its
    order (3548), in the quartic algebra of f = x^4-1000x^3-1000x^2-1000x-1000: exactly the ideal that PicardGroup
    tests for principality to glue Pic(O) to the Picard group of the equation order. Its Z-basis has entries with
    thousands of digits; SmallRepresentative multiplies by a short element of (O:cI^n) and LLL-reduces the HNF of the
    result, which has small entries even if a*ZBasis(cI^n) does not.
    Inputs (fixed, see README): cI as computed by CoprimeRepresentative with a2f3804 on 2026-11-10, stored as a
    Z-basis; the input of SmallRepresentative is the HNF basis of cI^3548, the same in every version.
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
// fixed inputs (the Z-basis of cI), printed with PrintSeqAlgEtQElt
data_cI:=[
    [<[1,-1,1/2,37/100]>,<[0,-1,-3/5,109/100]>,<[-1,0,1,9/20]>,<[-1,-1,-9/10,-89/50]>]
];
f:=x^4-1000*x^3-1000*x^2-1000*x-1000;
A:=EtaleAlgebra(f);
O:=MaximalOrder(A);
GO,gO:=PicardGroup(O : GRH:=true);
cI:=Ideal(O,[ A!t : t in data_cI[1] ]);
assert Order(cI@@gO) eq 3548;
II:=cI^3548;
zbII:=ZBasis(Ideal(O,ZBasis(II))); // the HNF basis, the same in every version
chk:=checksum(zbII);
printf "### INPUT %o\n",chk;
assert chk eq 199545024928899935;
t:=0;
for i in [1..5] do
    J:=Ideal(O,zbII); // a fresh copy
    J`ZBasis:=zbII;
    t0:=Cputime();
        aJ,a:=SmallRepresentative(J);
    t+:=Cputime(t0);
    assert aJ eq a*II;
    assert aJ subset OneIdeal(O);
end for;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
