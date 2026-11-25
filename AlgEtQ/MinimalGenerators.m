/////////////////////////////////////////////////////
// Copyright 2025.
// Stefano Marseglia, stefano.marseglia89@gmail.com
// https://stmar89.github.io/index.html
//
// Distributed under the terms of the CC-BY 4.0 licence.
// https://creativecommons.org/licenses/by/4.0/
/////////////////////////////////////////////////////



freeze;

//------------
// TwoGeneratingSet
//------------
import "ShortEltSmallRep.m" : short_non_zero_divisors;

///## Two-generating set

/// A procedure that given an invertible fractional ideal $I$, stores in the attribute `Generators` of $I$ one or two non-zerodivisors in $I$ generating $I$. The first one is `ShortElement(I)`, and the second one is chosen by a deterministic process: the output does not depend on the random state.
intrinsic TwoGeneratingSet(I::AlgEtQIdl)
{A procedure that given an invertible fractional ideal I, stores in the attribute Generators of I one or two non-zerodivisors in I generating I. The first one is ShortElement(I), and the second one is chosen by a deterministic process: the output does not depend on the random state.}
    require IsInvertible(I) : "the ideal must be invertible";
    if #Generators(I) gt 2 then
        a:=ShortElement(I);
        S:=Order(I);
        zbaS:=[a*z : z in ZBasis(S)];
        Q,q:=Quotient(I,zbaS);
        if IsTrivial(Q) then
            I`Generators:=[a]; //the ideal is principal
        else
            aS:=Ideal(S,zbaS);
            aS`ZBasis:=zbaS;
            // a and b generate I if and only if the images of b*ZBasis(S) generate Q
            generates:=func< b | Q eq sub<Q|[q(b*z):z in ZBasis(S)]> >;
            // We try elements of I in a fixed order: the LLL-reduced ZBasis of I (computed by
            // ShortElement) and small combinations of it, see short_non_zero_divisors. We take b
            // short modulo aS: the generators are used to multiply ideals, which is faster with
            // small generators.
            found:=false;
            for c in short_non_zero_divisors(ZBasis(I),30) do
                b:=ShortElementModulo(c,aS);
                if generates(b) then
                    found:=true;
                    break;
                end if;
            end for;
            if not found then
                // This can happen when Q has many local components with small residue fields. Then
                // b is built with CoprimeRepresentative: since I is invertible, b*S+a*S=I if and
                // only if b*I^-1+a*I^-1=S, that is, if and only if b*I^-1 is coprime to the
                // integral ideal a*I^-1. CoprimeRepresentative(I^-1,a*I^-1) returns such b, in
                // (I^-1)^-1=I.
                vprintf AlgEtQIdl,2: "TwoGeneratingSet: no candidate works, we use CoprimeRepresentative\n";
                Iinv:=Inverse(I);
                b:=ShortElementModulo(CoprimeRepresentative(Iinv,a*Iinv),aS);
                assert2 generates(b);
            end if;
            // b+t*a has the same class modulo aS; it is a zero-divisor for at most one t per
            // component
            t:=0;
            while IsZeroDivisor(b+t*a) do
                t+:=1;
            end while;
            b:=b+t*a;
            I`Generators:=[a,b];
        end if;
    end if;
end intrinsic;

/* TESTS

    printf "### Testing MinimalGenerators:";
    //AttachSpec("../spec");
    SetAssertions(2);
    SetClassGroupBounds("GRH");
    _<x>:=PolynomialRing(Integers());
    f:=x^4-1000*x^3-1000*x^2-1000*x-1000;
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    SetAssertions(1); // with SetAssertions(2), PicardGroup takes ~100 secs for this order: here we
                      // test MinimalGenerators
    P,p:=PicardGroup(E : GRH:=true);
    SetAssertions(2);

    for g in Generators(P) do
        I:=p(g);
        TwoGeneratingSet(I);
        assert #Generators(I) le 2;
        printf ".";
    end for;

    // The second generator is short modulo a*E, where a is the first one.
    for g in Generators(P) do
        I:=Ideal(E,ZBasis(p(g))); // with 4 generators, so that TwoGeneratingSet computes new ones
        TwoGeneratingSet(I);
        if #Generators(I) eq 2 then
            a,b:=Explode(Generators(I));
            zbaE:=[ a*z : z in ZBasis(E) ];
            aE:=Ideal(E,zbaE);
            aE`ZBasis:=zbaE;
            assert ShortElementModulo(b,aE) eq b;
            assert I eq Ideal(E,[a,b]);
        end if;
        printf ".";
    end for;


    // TwoGeneratingSet is deterministic: the same generators with another seed, also when no
    // candidate works and CoprimeRepresentative is used (in a product of four imaginary quadratic
    // fields where 2 splits, Q = I/aS has eight local components with residue field F_2 for I =
    // 2O*P).
    A2:=EtaleAlgebra([ NumberField(x^2-x+k) : k in [2,4,6,8] ]);
    O2:=MaximalOrder(A2);
    pp2:=PrimesAbove(2*O2);
    ids:=[ Ideal(E,ZBasis(p(g))) : g in Generators(P) ] cat [ Ideal(O2,ZBasis(2*Q)) : Q in pp2[1..3] ];
    for I in ids do
        S:=Order(I);
        I1:=Ideal(S,ZBasis(I));
        SetSeed(1);
        TwoGeneratingSet(I1);
        I2:=Ideal(S,ZBasis(I));
        SetSeed(2);
        TwoGeneratingSet(I2);
        assert Generators(I1) eq Generators(I2);
        assert I eq Ideal(S,Generators(I1)) and forall{ g : g in Generators(I1) | not IsZeroDivisor(g) };
        printf ".";
    end for;
    SetSeed(1);

    // test if TwoGeneratingSet makes the power faster
    // Conclusion: yes. By quite a bit!
    f:=x^4-100*x^3-100*x^2-100*x-100;
    A:=EtaleAlgebra(f);
    E:=EquationOrder(A);
    P,p:=PicardGroup(E : GRH:=true);
    repeat
        Ii:=Random(P);
    until Ii ne Zero(P);
    I:=p(Ii);

    delete I`IsInvertible;
    exp:=[ Random(2,30) : i in [1..100]];
    l1:=[ I^i : i in exp ];
    printf ".";

    assert IsInvertible(I);
    TwoGeneratingSet(I);
    assert #Generators(I) eq 2;
    l2:=[ I^i : i in exp ];
    assert l1 eq l2;
    printf ".";

    I:=SmallRepresentative(I);
    delete I`IsInvertible;
    l1:=[ I^i : i in exp ];

    assert IsInvertible(I);
    TwoGeneratingSet(I);
    l2:=[ I^i : i in exp ];
    assert l1 eq l2;
    SetAssertions(1);
    printf " all good!";
*/