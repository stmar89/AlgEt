/////////////////////////////////////////////////////
// Copyright 2025.
// Stefano Marseglia, stefano.marseglia89@gmail.com
// https://stmar89.github.io/index.html
//
// Distributed under the terms of the CC-BY 4.0 licence.
// https://creativecommons.org/licenses/by/4.0/
/////////////////////////////////////////////////////



freeze;

declare attributes AlgEtQIdl : InertiaDegree,
                               RamificationIndex;

declare verbose Uniformizers, 3;

// we continue the doc of FactPrimes.m

///## Valuations, and primes of the maximal order

valuation_at_invertible_prime:=function(J,P)
// Given an integral ideal J and an invertible prime P of the same order S, returns the valuation of
// J at P. The localization S_P is a discrete valuation ring, so J_P = P^v S_P. For k >= v, J+P^k is
// P^v at P and S at every other prime Q (P^k is not contained in Q), so [S : J+P^k] = p^(f*v),
// where [S:P] = p^f. An upper bound for v: [S:J] >= [S_P:J_P] = p^(f*v).
    if not J subset P then
        return 0;
    end if;
    S:=Order(P);
    p:=MinimalInteger(P); // P meet Z = pZ
    f:=Ilog(p,Integers()!Index(S,P));
    k:=Valuation(Integers()!Index(S,J),p) div f;
    v:=Valuation(Integers()!Index(S,J+P^k),p) div f;
    assert2 (J+P^(k+1)) eq (J+P^k); // v < k+1: the bound holds
    return v;
end function;

/// Given an element $x$ and a prime $P$ of an order $S$, returns the valuation of $x$ at $P$. Either $S$ is the maximal order, or $P$ is invertible and $x$ is not a zero-divisor; in the second case the maximal order is not computed.
intrinsic Valuation(x::AlgEtQElt,P::AlgEtQIdl) -> RngIntElt
{Given an element x and a prime P of an order S, returns the valuation of x at P. Either S is the maximal order, or P is invertible and x is not a zero-divisor; in the second case the maximal order is not computed.}
    OO:=Order(P);
    require IsPrime(P): "P must be a prime";
    if (assigned OO`IsMaximal and OO`IsMaximal) or (assigned Algebra(OO)`MaximalOrder and IsMaximal(OO)) then
        _,comp:=IsProductOfIdeals(P);
        ind:=[ i : i in [1..#comp] | not One(NumberField(Order(comp[i]))) in comp[i] ];
        assert #ind eq 1;
        ind:=ind[1];
        return Valuation(Components(x)[ind],comp[ind]);
    end if;
    require IsInvertible(P): "P must be a prime of the maximal order, or an invertible prime";
    require not IsZeroDivisor(x): "x must not be a zero-divisor, unless P is a prime of the maximal order";
    return Valuation(x*OneIdeal(OO),P);
end intrinsic;

/// Given a fractional ideal $I$ and a prime $P$ of the same order $S$, returns the valuation of $I$ at $P$. Either $S$ is the maximal order, or $P$ is invertible, so that the localization of $S$ at $P$ is a discrete valuation ring; in the second case the maximal order is not computed.
intrinsic Valuation(I::AlgEtQIdl,P::AlgEtQIdl) -> RngIntElt
{Given a fractional ideal I and a prime P of the same order S, returns the valuation of I at P. Either S is the maximal order, or P is invertible, so that the localization of S at P is a discrete valuation ring; in the second case the maximal order is not computed.}
    OO:=Order(P);
    require Order(I) eq OO and IsPrime(P): "P must be a prime of the order of I";
    if (assigned OO`IsMaximal and OO`IsMaximal) or (assigned Algebra(OO)`MaximalOrder and IsMaximal(OO)) then
        _,compI:=IsProductOfIdeals(I);
        _,compP:=IsProductOfIdeals(P);
        ind:=[ i : i in [1..#compP] | not One(NumberField(Order(compP[i]))) in compP[i] ];
        assert #ind eq 1;
        ind:=ind[1];
        return Valuation(compI[ind],compP[ind]);
    end if;
    require IsInvertible(P): "P must be a prime of the maximal order, or an invertible prime";
    // d*I is integral: v_P(I) = v_P(d*I) - v_P(d), with v_P(d) = v_p(d)*v_P(p)
    J,d:=MakeIntegral(I);
    v:=valuation_at_invertible_prime(J,P);
    p:=MinimalInteger(P);
    if d mod p eq 0 then
        v -:= Valuation(d,p)*valuation_at_invertible_prime(p*OneIdeal(OO),P);
    end if;
    return v;
end intrinsic;

/// Given a prime $P$ of the maximal order $\mathcal{O}$ above the rational prime $p$, it returns the inertia degree of $P$, that is, the index of the finite field extension $\mathbb{F}_p \to \mathcal{O}/P$.
intrinsic InertiaDegree(P::AlgEtQIdl) -> RngIntElt
{Given a prime P of the maximal order O above the rational prime p, it returns the inertia degree of P, that is, the index of the finite field extension GF(p)->O/P.}
    if not assigned P`InertiaDegree then
        OO:=Order(P);
        require IsMaximal(OO) and IsPrime(P): "The input needs to be a prime of the maximal order.";
        q:=Integers()!Index(OO,P);
        p:=MinimalInteger(P); // P meet Z = pZ
        f:=Ilog(p,q);
        assert p^f eq q;
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
        p:=MinimalInteger(P); // P meet Z = pZ
        P`RamificationIndex:=Valuation(Algebra(P)!p,P);
    end if;
    return P`RamificationIndex;
end intrinsic;

/// Given a sequence of primes $P$ of the maximal order, it returns a sequence of elements $t_P$ such that $t_P$ is a uniformizer at $P$ and a unit at every other prime in the sequence. The output does not depend on the random state.
intrinsic Uniformizers(PPs::SeqEnum[AlgEtQIdl]) -> SeqEnum
{Given a sequence of primes P of the maximal order, it returns a sequence of elements t_P such that t_P is a uniformizer at P and a unit at every other prime in the sequence. The output does not depend on the random state.}
    require #PPs gt 0 : "The sequence must not be empty.";
    OO:=Order(PPs[1]);
    require IsMaximal(OO) and forall{P : P in PPs | Order(P) eq OO and IsPrime(P)} : "The input needs to be a sequence of primes of the maximal order.";
    require #Seqset(PPs) eq #PPs : "The primes must be distinct.";
    vprintf Uniformizers,2 : "nice_uniformizers: %o primes\n",#PPs;
    one:=One(Algebra(PPs[1]));
    nice_unifs:=[];
    PPs2:=[ P^2 : P in PPs];
    // Only the classes of t_P modulo the P_i^2 matter, so t_P is replaced at the end by a short
    // element of its class modulo L, the product of the P_i^2 (their intersection, since they are
    // coprime): the LLL reduction of the HNF of L, so that the output depends only on the primes.
    L:=&meet PPs2;
    mL:=MinimalInteger(L);
    for iP->P in PPs do
        vprintf Uniformizers,2 : "nice_uniformizers: %o-th prime...",iP;
        P2:=PPs2[iP];
        // An element of P not in P^2 has valuation 1 at P. We take the first one in the LLL-reduced
        // Z-basis of P obtained from its HNF, which depends only on P: some element of a Z-basis of
        // P is not in P^2, since P ne P^2.
        tP:=[ z : z in ZBasisLLLFromHNF(P) | not z in P2 ][1];
        vprintf Uniformizers,2 : "got an element of valuation 1...";
        assert2 tP in P and not tP in P2;
        elts:=[ i eq iP select tP else one : i in [1..#PPs] ];
        tP:=CRT(PPs2,elts);
        tP:=ShortElementModulo(tP,L : Canonical:=true);
        // Adding a multiple of the smallest positive integer mL in L does not change the classes
        // modulo the P_i^2, and tP+k*mL is a zero-divisor for at most one k per component.
        k:=0;
        while IsZeroDivisor(tP+k*mL) do
            k+:=1;
        end while;
        tP:=tP+k*mL;
        vprintf Uniformizers,2 : "got a short non-zerodivisor...";
        vprintf Uniformizers,2 : "asserts";
        assert2 tP in P and not tP in P2;
        vprintf Uniformizers,2 : ".";
        assert2 tP in MaximalOrder(Algebra(tP));
        vprintf Uniformizers,2 : ".";
        assert2 forall{ Q : Q in [ Q : Q in PPs | Q ne P ] | not tP in Q };
        vprintf Uniformizers,2 : ".";
        // for some reason this is very time consuming
        assert3 Seqset(PrimesAbove(tP*MaximalOrder(Algebra(tP)))) meet Seqset(PPs) eq {P};
        vprintf Uniformizers,2 : ".";
        assert not IsZeroDivisor(tP);
        vprintf Uniformizers,2 : "all passed\n";
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
            // short: the reduction modulo L, the product of the P^2, plus k times the smallest
            // integer of L, with k at most the number of components
            L:=&meet[ P^2 : P in pp ];
            assert forall{ t : t in ts | exists{ k : k in [0..#Components(A)] |
                                            t-ShortElementModulo(t,L : Canonical:=true) eq k*MinimalInteger(L) } };
            assert forall{ 0 : i,j in [1..#pp] | Valuation(ts[i],pp[j]) eq (i eq j select 1 else 0) };
            // not random: the same uniformizers with another seed, for equal primes created again
            SetSeed(2);
            assert Uniformizers([ Ideal(O,ZBasis(P)) : P in pp ]) eq ts;
            SetSeed(1);
        end for;
        printf ".";
    end for;
    // the primes must be distinct
    _<x>:=PolynomialRing(Integers());
    A:=EtaleAlgebra(x^2+1); O:=MaximalOrder(A);
    P:=PrimesAbove(2*OneIdeal(O))[1];
    failed:=false; try _:=Uniformizers([P,P]); catch e failed:=true; end try; assert failed;
    printf ".";
    // the primes must be of the same maximal order
    _<x>:=PolynomialRing(Integers());
    O1:=MaximalOrder(EtaleAlgebra(x^2+1)); O2:=MaximalOrder(EtaleAlgebra(x^2+1));
    failed:=false; try _:=Uniformizers([PrimesAbove(2*OneIdeal(O1))[1],PrimesAbove(5*OneIdeal(O2))[1]]); catch e failed:=true; end try; assert failed;
    failed:=false; try _:=Uniformizers([Universe([PrimesAbove(2*OneIdeal(O1))[1]])|]); catch e failed:=true; end try; assert failed;
    printf ".";
    // Valuation at an invertible prime P of a non-maximal order, without the maximal order: the
    // same as the valuation at the prime P*O of the maximal order, for fractional ideals (also not
    // coprime with the conductor) and elements; a singular prime is refused
    _<x>:=PolynomialRing(Integers());
    for f in [ (x^2+5)*(x^2+7), x^6+3*x^5+99*x^2+2997 ] do
        A:=EtaleAlgebra(f); E:=EquationOrder(A); a:=PrimitiveElement(A);
        pp:=&cat[ PrimesAbove(q*OneIdeal(E)) : q in [2,3,5,7] ];
        reg:=[ P : P in pp | IsInvertible(P) ];
        sing:=[ P : P in pp | not IsInvertible(P) ];
        ids:=[ (1/6)*reg[1]^2*reg[#reg], reg[2]^3*(2*a+1)*OneIdeal(E), (1/35)*OneIdeal(E) ] cat [ P*Q : P in sing, Q in reg[1..2] ];
        elts:=[ a, 2*a+1, (a^2-3)/10, A!12 ];
        vals:=[ [ Valuation(I,P) : P in reg ] : I in ids ];
        evals:=[ [ Valuation(y,P) : P in reg ] : y in elts ];
        assert not assigned A`MaximalOrder;
        O:=MaximalOrder(A);
        regO:=[ O!!P : P in reg ]; // the primes of O above the regular primes
        assert forall{ P : P in regO | IsPrime(P) };
        assert vals eq [ [ Valuation(O!!I,P) : P in regO ] : I in ids ];
        assert evals eq [ [ Valuation(y,P) : P in regO ] : y in elts ];
        if #sing gt 0 then
            try _:=Valuation(ids[1],sing[1]); ok:=false; catch e ok:=true; end try;
            assert ok;
        end if;
    end for;
    printf ".";
    SetAssertions(1);
    printf " all good!";

*/
