/////////////////////////////////////////////////////
// Copyright 2025.
// Stefano Marseglia, stefano.marseglia89@gmail.com
// https://stmar89.github.io/index.html
//
// Distributed under the terms of the CC-BY 4.0 licence.
// https://creativecommons.org/licenses/by/4.0/
/////////////////////////////////////////////////////



freeze;

declare verbose ShortEltSmallRep, 2;

declare attributes AlgEtQIdl : ShortElement, SmallRepresentative;

///# Short Elements
/// Let $I$ be a fractional ideal in an étale algebra $A$ over $\mathbb{Q}$.
/// By LLL-reducing a $\mathbb{Z}$-basis of $I$, written in the coordinates of the absolute basis of $A$, and then picking an element with small coefficients with respect to that basis allows producing elements with small coefficients with respect to the absolute basis of $A$.

//------------
// ShortElement
//------------

short_non_zero_divisors:=function(zb,N)
// Given a short Z-basis zb of a fractional ideal, returns short elements of the ideal which are not
// zero-divisors, in a fixed order: the elements of zb which are not zero-divisors, then the ones
// among N combinations of the elements of zb with coefficients in [-3..3] from a Lehmer generator
// with a fixed initial state. The combinations are needed in products of number fields, where the
// elements of zb can all be zero-divisors (for example for an ideal of the maximal order). The
// output depends only on zb: it is not random.
    cands:=[ z : z in zb | not IsZeroDivisor(z) ];
    s:=1;
    for k in [1..N] do
        cs:=[];
        for i in [1..#zb] do
            s:=(48271*s) mod (2^31-1);
            Append(~cs,(s mod 7)-3);
        end for;
        z:=DotProduct(cs,zb);
        if not IsZeroDivisor(z) then
            Append(~cands,z);
        end if;
    end for;
    return cands;
end function;

/// Given a fractional ideal $I$, returns a non-zerodivisor in $I$ with small coefficients (in the LLL sense). This is achieved by picking, among the non-zerodivisors in an LLL-reduced $\mathbb{Z}$-basis of $I$ (in the absolute coordinates), one whose norm has the smallest absolute value. If all the elements of this basis are zerodivisors, then combinations of them with small coefficients, in a fixed order, are used instead: the output is not random. If the stored ZBasis is not LLL-reduced, it will be changed to an LLL-reduced one.
intrinsic ShortElement(I::AlgEtQIdl) ->AlgEtQElt
{Given a fractional ideal I, returns a non-zerodivisor in I with small coefficients (in the LLL sense). This is achieved by picking, among the non-zerodivisors in an LLL-reduced Z-basis of I (in the absolute coordinates), one whose norm has the smallest absolute value. If all the elements of this basis are zerodivisors, then combinations of them with small coefficients, in a fixed order, are used instead: the output is not random. If the stored ZBasis is not LLL-reduced, it will be changed to an LLL-reduced one.}
    if not assigned I`ShortElement then
        ZBasisLLL(I);
        B := ZBasis(I);
        elts:=[ x : x in B | not IsZeroDivisor(x) ];
        // If all the elements of B are zero-divisors, which happens in products of number fields
        // (for example for ideals of the maximal order), we use combinations of them (see
        // short_non_zero_divisors): 30 are almost always enough, otherwise we take more.
        N:=30;
        while #elts eq 0 do
            elts:=short_non_zero_divisors(B,N);
            N*:=2;
        end while;
        _,i:=Min([Abs(Norm(x)) : x in elts]);
        elt:=elts[i];
        I`ShortElement:= elt;
    end if;
    return I`ShortElement;
end intrinsic;

//------------
// SmallRepresentative
//------------

/// Given a fractional $R$-ideal $I$, returns the fractional $R$-ideal $a*I$, and the element $a$, where $a$ is a non-zero divisor chosen such that $a\cdot I \subseteq  R$, and the cardinality of $R/aI$ is small. More precisely, $a$ is the output of `ShortElement` when applied to $(R:I)$. Note that if $I$ is invertible then $R/aI$ is isomorphic to $(R:I)/aR$.
intrinsic SmallRepresentative(I::AlgEtQIdl) ->AlgEtQIdl,AlgEtQElt
{Given a fractional R-ideal I, returns the fractional R-ideal a*I, and the element a, where a is chosen so that a*I is a subset of R, and the cardinality of R/aI is small. More precisely, a is the output of ShortElement when applied to (R:I). Note that if I is invertible R/aI is isomorphic to (R:I)/aR.}
    if not assigned I`SmallRepresentative then
        R:=Order(I);
        // if I is invertible and its inverse is known, then the inverse is (R:I)
        cRI:=assigned I`Inverse select I`Inverse else ColonIdeal(R,I);
        a:=ShortElement(cRI);
        aI:=a*I;
        // The ZBasis of aI might be very big. We make it smaller.
        // Since aI is contained in R with small index, its HNF has small entries even if ZBasis(I)
        // is huge, and ZBasisLLL will use it: this is much faster than the LLL of a*ZBasis(I).
        ZBasisLLL(aI);
        // If I is known to be invertible, then the inverse of aI is free: (aI)^-1 = a^-1*I^-1, and
        // I^-1 = (R:I). The two ideals are stored as inverses of each other, as Inverse does.
        if assigned I`IsInvertible and I`IsInvertible then
            aI`IsInvertible:=true;
            aI`MultiplicatorRing:=R;
            aIinv:=(1/a)*cRI;
            aIinv`IsInvertible:=true;
            aIinv`MultiplicatorRing:=R;
            aIinv`Inverse:=aI;
            aI`Inverse:=aIinv;
        end if;
        vprintf ShortEltSmallRep,2: "SmallRepresentative:\n
                                I = %o\n,aI = %o\n",PrintSeqAlgEtQElt(ZBasis(I)),PrintSeqAlgEtQElt(ZBasis(aI));
        I`SmallRepresentative:=<aI,a>;
    end if;
    return Explode(I`SmallRepresentative);
end intrinsic;


/* TESTS

    printf "### Testing ShortEltSmallRep:";
    //AttachSpec("../spec");
    SetAssertions(2);
    _<x>:=PolynomialRing(Integers());
    f:=(x^2+5)*(x^2+7)*(x^2+11);
    assert IsSquarefree(f);
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    ff:=Conductor(E);
    _:=ShortElement(ff);
    oo:=FindOverOrders(E);
    for S in oo do
        printf ".";
        ff:=Conductor(S);
        _:=ShortElement(ff);
    end for;
    // SmallRepresentative: aI = a*I is integral, a is a non-zero-divisor in (R:I), the output is
    // stored, and for invertible I the index of aI in R is the index of a*R in (R:I) = I^-1.
    for R in [ E ] cat oo[1..Min(5,#oo)] do
        cR:=Conductor(R);
        for I in [ cR, OneIdeal(R), 7*cR, (1/3)*cR, cR^2 ] do
            aI,a:=SmallRepresentative(I);
            assert aI eq a*I and IsIntegral(aI) and not IsZeroDivisor(a) and a in ColonIdeal(R,I);
            aI2,a2:=SmallRepresentative(I);
            assert IsIdentical(aI2,aI) and a2 eq a;
            if IsInvertible(I) then
                assert Index(OneIdeal(R),aI) eq Index(ColonIdeal(R,I),a*OneIdeal(R));
            end if;
        end for;
        printf ".";
    end for;

    // SmallRepresentative of an ideal known to be invertible: the inverse of aI is stored, it is
    // (R:aI), and its stored inverse is aI.
    s:=1;
    for R in [ E ] cat oo[1..Min(5,#oo)] do
        AR:=Algebra(R);
        for k in [1..4] do
            gs:=[];
            for j in [1..2] do
                cs:=[];
                for i in [1..AbsoluteDimension(AR)] do
                    s:=(48271*s) mod (2^31-1);
                    Append(~cs,(s mod 7)-3);
                end for;
                Append(~gs,DotProduct(cs,ZBasis(R)));
            end for;
            if IsZeroDivisor(gs[1]) then
                continue;
            end if;
            I:=Ideal(R,gs);
            if IsInvertible(I) then
                aI,a:=SmallRepresentative(I);
                assert assigned aI`Inverse and aI`Inverse eq ColonIdeal(R,aI);
                assert aI*aI`Inverse eq OneIdeal(R) and IsIdentical(aI`Inverse`Inverse,aI);
            end if;
        end for;
        printf ".";
    end for;

    // ShortElement is not random: in a product of number fields, the LLL-reduced Z-basis of an
    // ideal of the maximal order consists of zero-divisors, and combinations of it are used. The
    // output is the same for another seed.
    A:=EtaleAlgebra([ NumberField(x^2+1), NumberField(x^2-2), NumberField(x^3-5) ]);
    O:=MaximalOrder(A);
    for I in [ OneIdeal(O), 3*OneIdeal(O), &*PrimesAbove(6*O) ] do
        ZBasisLLL(I);
        SetSeed(1);
        I1:=Ideal(O,ZBasis(I));
        a:=ShortElement(I1);
        assert a in I and not IsZeroDivisor(a);
        SetSeed(2);
        I2:=Ideal(O,ZBasis(I));
        assert ShortElement(I2) eq a;
        printf ".";
    end for;
    SetSeed(1);
    SetAssertions(1);
    printf " all good!";
*/
