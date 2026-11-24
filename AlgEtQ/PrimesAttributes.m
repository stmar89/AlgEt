/////////////////////////////////////////////////////
// Copyright 2025.
// Stefano Marseglia, stefano.marseglia89@gmail.com
// https://stmar89.github.io/index.html
//
// Distributed under the terms of the CC-BY 4.0 licence.
// https://creativecommons.org/licenses/by/4.0/
/////////////////////////////////////////////////////



declare attributes AlgEtQIdl : InertiaDegree,
                               RamificationIndex;

declare verbose nice_unif,3;

// we continue the doc of FactPrimes.m

///## Primes of the maximal order

/// Given an element $x$ and a prime $P$ of the maximal order, returns the valuation of $x$ at $P$.
intrinsic Valuation(x::AlgEtQElt,P::AlgEtQIdl) -> RngIntElt
{Given an element x and a prime P of the maximal order, returns the valuation of x at P.}
    OO:=Order(P);
    require IsMaximal(OO) and IsPrime(P): "The input needs to be a prime of the maximal order.";
    _,comp:=IsProductOfIdeals(P);
    ind:=[ i : i in [1..#comp] | not One(NumberField(Order(comp[i]))) in comp[i] ];
    assert #ind eq 1;
    ind:=ind[1];
    return Valuation(Components(x)[ind],comp[ind]);
end intrinsic;

/// Given a fractional ideal $I$ and a prime $P$, both of the maximal order, returns the valuation of $I$ at $P$.
intrinsic Valuation(I::AlgEtQIdl,P::AlgEtQIdl) -> RngIntElt
{Given a fractional ideal I and a prime P, both of the maximal order, returns the valuation of I at P.}
    OO:=Order(P);
    require IsMaximal(OO) and IsPrime(P): "The input needs to be a prime of the maximal order.";
    _,compI:=IsProductOfIdeals(I);
    _,compP:=IsProductOfIdeals(P);
    ind:=[ i : i in [1..#compP] | not One(NumberField(Order(compP[i]))) in compP[i] ];
    assert #ind eq 1;
    ind:=ind[1];
    return Valuation(compI[ind],compP[ind]);
end intrinsic;

/// Given a prime $P$ of the maximal order $\mathcal{O}$ above the rational prime $p$, it returns the inertia degree of $P$, that is, the index of the finite field extension $\mathbb{F}_p \to \mathcal{O}/P$.
intrinsic InertiaDegree(P::AlgEtQIdl) -> RngIntElt
{Given a prime P of the maximal order O above the rational prime p, it returns the inertia degree of P, that is, the index of the finite field extension GF(p)->O/P.}
    if not assigned P`InertiaDegree then
        OO:=Order(P);
        require IsMaximal(OO) and IsPrime(P): "The input needs to be a prime of the maximal order.";
        q:=Index(OO,P);
        t,p,f:=IsPrimePower(q);
        assert t;
        P`InertiaDegree:= f;
    end if;
    return P`InertiaDegree;
end intrinsic;

/// Given a prime $P$ of the maximal order $\mathcal{O}$, it returns the ramification index of $P$.
intrinsic RamificationIndex(P::AlgEtQIdl) -> RngIntElt
{Given a prime P of the maximal order O, it returns the ramification index of P.}
    if not assigned P`RamificationIndex then
        OO:=Order(P);
        require IsMaximal(OO) and IsPrime(P): "The input needs to be a prime of the maximal order.";
        q:=Index(OO,P);
        t,p:=IsPrimePower(q);
        assert t;
        P`RamificationIndex:=Valuation(Algebra(P)!p,P);
    end if;
    return P`RamificationIndex;
end intrinsic;

/// Given a sequence of primes $P$ of the maximal order, it returns a sequence of elements $t_P$ such that $t_P$ is a uniformizer at $P$ and a unit at every other prime in the sequence.
intrinsic Uniformizers(PPs::SeqEnum[AlgEtQIdl]) -> SeqEnum
{Given a sequence of primess P of the maximal order, it returns a sequence of elements t_P such that t_P is a uniformizer at P and a unit at every other prime in the sequence.}
    OO:=Order(PPs[1]);
    require IsMaximal(OO) and forall{P : P in PPs | IsPrime(P)} : "The input needs to be a sequence of primess of the maximal order.";
    vprintf nice_unif,2 : "nice_uniformizers: %o primes\n",#PPs;
    one:=One(Algebra(PPs[1]));
    nice_unifs:=[];
    PPs2:=[ P^2 : P in PPs];
    for iP->P in PPs do
        vprintf nice_unif,2 : "nice_uniformizers: %o-th prime...",iP;
        P2:=PPs2[iP];
        Q,q:=Quotient(P,P2);
        repeat
            t:=Random(Q);
        until t ne Zero(Q);
        vprintf nice_unif,2 : "got non-zero elt...";
        tP:=t@@q;
        while IsZeroDivisor(tP) do
            tP +:=Random(P2);
        end while;
        vprintf nice_unif,2 : "got non-zerodivisor...";
        assert2 tP in P and not tP in P2;
        elts:=[ i eq iP select tP else one : i in [1..#PPs] ];
        tP:=CRT(PPs2,elts);
        vprintf nice_unif,2 : "asserts";
        assert2 tP in P and not tP in P2;
        vprintf nice_unif,2 : ".";
        assert2 tP in MaximalOrder(Algebra(tP));
        vprintf nice_unif,2 : ".";
        assert2 forall{ Q : Q in [ Q : Q in PPs | Q ne P ] | not tP in Q };
        vprintf nice_unif,2 : ".";
        // for some reason this is very time consuming
        assert3 Seqset(PrimesAbove(tP*MaximalOrder(Algebra(tP)))) meet Seqset(PPs) eq {P};
        vprintf nice_unif,2 : ".";
        assert not IsZeroDivisor(tP);
        vprintf nice_unif,2 : "all passed\n";
        Append(~nice_unifs,tP);
    end for;
    return nice_unifs;
end intrinsic;

/* TESTS

    printf "### Testing PrimesAttributes:";
    //AttachSpec("../spec");
    SetAssertions(2);
    _<x>:=PolynomialRing(Integers());
    // For the primes P of the maximal order above p: sum of e_P*f_P is the dimension, [O:P] =
    // p^f_P, the valuation of p at P is e_P, valuations are additive, and the uniformizers have
    // valuation 1 at their prime and 0 at the others.
    for f in [ x^4-1000*x^3-1000*x^2-1000*x-1000, (x^2+5)*(x^3-2000), (x-1)*(x^2+7) ] do
        A:=EtaleAlgebra(f);
        O:=MaximalOrder(A);
        F:=PrimitiveElement(A);
        elts:=[ F+1, F^2-3, 2*F+5 ];
        elts:=[ y : y in elts | not IsZeroDivisor(y) ];
        for p in [2,3,5,7] do
            pp:=PrimesAbove(p*O);
            assert &+[ RamificationIndex(P)*InertiaDegree(P) : P in pp ] eq AbsoluteDimension(A);
            for P in pp do
                assert Index(O,P) eq p^InertiaDegree(P);
                assert Valuation(A!p,P) eq RamificationIndex(P);
                assert Valuation(p*O,P) eq RamificationIndex(P);
                assert forall{ 0 : y,z in elts | Valuation(y*z,P) eq Valuation(y,P)+Valuation(z,P) };
                assert forall{ y : y in elts | Valuation(y*OneIdeal(O),P) eq Valuation(y,P) };
            end for;
            ts:=Uniformizers(pp);
            assert forall{ 0 : i,j in [1..#pp] | Valuation(ts[i],pp[j]) eq (i eq j select 1 else 0) };
        end for;
        printf ".";
    end for;
    SetAssertions(1);
    printf " all good!";

*/
