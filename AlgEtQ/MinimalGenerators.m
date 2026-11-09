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
import "Ord.m" : MatrixAtoQ,MatrixAtoZ;

///## Two-generating set

/// A procedure that given an invertible fractional ideal $I$, stores in the attribute `Generators` of $I$ one or two non-zerodivisors in $I$ generating $I$. The elements are chosen by a non-deterministic process.
intrinsic TwoGeneratingSet(I::AlgEtQIdl)
{A procedure that given an invertible fractional ideal I, stores in the attribute Generators of I one or two non-zerodivisors in I generating I. The elements are chosen by a non-deterministic process.}
    require IsInvertible(I) : "the ideal must be invertible";
    if #Generators(I) gt 2 then
        // if assigned I`IsIntegral and IsIntegral(I) then
        //     a:=MinimalInteger(I);
        // else
        //     a:=Random(I : ZeroDivisorsAllowed:=false );
        // end if;
        a:=ShortElement(I);
        S:=Order(I);
        zbaS:=[a*z : z in ZBasis(S)];
        Q,q:=Quotient(I,zbaS);
        if IsTrivial(Q) then
            I`Generators:=[a]; //the ideal is principal
        else
            aS:=Ideal(S,zbaS);
            aS`ZBasis:=zbaS;
            repeat
                repeat
                    b:=Random(Q);
                until b ne Zero(Q);
                // We take b short modulo aS: the generators are used to multiply ideals, which is faster with small generators.
                b:=ShortElementModulo(b@@q,aS);
            //until I eq Ideal(S,[a,b]);
            until Q eq sub<Q|[q(b*z):z in ZBasis(S)]>;
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
    SetAssertions(1); // with SetAssertions(2), PicardGroup takes ~100 secs for this order: here we test MinimalGenerators
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