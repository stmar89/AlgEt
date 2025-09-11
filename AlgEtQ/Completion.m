/////////////////////////////////////////////////////
// Copyright 2025.
// Stefano Marseglia, stefano.marseglia89@gmail.com
// https://stmar89.github.io/index.html
//
// Distributed under the terms of the CC-BY 4.0 licence.
// https://creativecommons.org/licenses/by/4.0/
/////////////////////////////////////////////////////



freeze;

/// Given a prime $P$ of the maximal order of an étale algebra $A$, returns the $p$-adic field corresponding to the completion $A_P$ and the natural homomorphism $A\to A_P$ (with preimages). The parameter `MinPrecision` is passed to `Completion`.
intrinsic Completion(P::AlgEtQIdl : MinPrecision:=20) -> FldPad,Map
{Given a prime of the maximal order of an etale algebra A, returns the p-adic field corresponding to the completion A_P and the natural homormophism A->A_P (with preimages). The parameter MinPrecision is passed to Completion.}
    L:=Algebra(P);
    require IsMaximal(Order(P)) and IsPrime(P) : "the ideal must be a prime ideal of the maximal order";
    nfs,embs:=Components(L);
    test,PKs:=IsProductOfIdeals(P);
    assert test;
    ind:=[ i : i in [1..#PKs] | not (Order(PKs[i]) ! 1) in PKs[i] ];
    assert #ind eq 1;
    ind:=ind[1];
    K:=nfs[ind];
    mK:=embs[ind]; // mK:=K->L
    PK:=PKs[ind];
    LP,mLP:=Completion(K,PK : Precision:=MinPrecision); // mLP:K->LP
    // it seems that Magma Ignore the Parameter Precision and returns always a ring with Infinity precition.
    // workaround
    LP1:=ChangePrecision(LP,MinPrecision);
    map:=hom< L->LP1 | x:->LP1!mLP(Components(x)[ind]) ,
                      y:->mK((LP!y)@@mLP) >;
    return LP1,map;
end intrinsic;

///# Example 4
/// ```
/// _<x>:=PolynomialRing(Integers());
/// f:=(x^8+16)*(x^8+81);
/// A:=EtaleAlgebra(f);
/// O:=MaximalOrder(A);
/// // Consider a bunch of prime of O and their uniformizers in O.
/// pp:=PrimesAbove(2*3*5*7*O);
/// unifs:=Uniformizers(pp);
/// // We now verify that each element is a uniformizer at the correct prime and a unit everywhere else
/// for iP->P in pp do
///     AP,mP:=Completion(P);
///     iP,[ Valuation(mP(t)) : t in unifs ];
/// end for;
/// ```


/* TESTS

    printf "### Testing Completion:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    PP<x>:=PolynomialRing(Integers());
    polys:=[
        x^6+3*x^4-10*x^3+15*x^2+125,
        (x^2+5)*(x^4-4*x^3+5*x^2-20*x+25),
        (x^4-5*x^3+15*x^2-25*x+25)*(x^4+5*x^3+15*x^2+25*x+25)
        ];
    for h in polys do
        L:=EtaleAlgebra(h);
        a:=PrimitiveElement(L);
        O:=MaximalOrder(L);
        p:=5;
        pp:=PrimesAbove(p*O);
        for P in pp do
            C,mC:=Completion(P);
        end for;
        printf ".";
    end for;

    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    O:=MaximalOrder(A);
    // Consider a bunch of prime of O and their uniformizers in O.
    pp:=PrimesAbove(2*3*5*7*O);
    unifs:=Uniformizers(pp);
    // We now verify that each element is a uniformizer at the correct prime and a unit everywhere else
    for iP->P in pp do
        AP,mP:=Completion(P);
        for it->t in unifs do
            if iP eq it then
                assert Valuation(mP(t)) eq 1;
            else
                assert Valuation(mP(t)) eq 0;
            end if;
        end for;
    end for;

    printf " all good!";
*/