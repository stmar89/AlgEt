/////////////////////////////////////////////////////
// Copyright 2025.
// Stefano Marseglia, stefano.marseglia89@gmail.com
// https://stmar89.github.io/index.html
// 
// Distributed under the terms of the CC-BY 4.0 licence.
// https://creativecommons.org/licenses/by/4.0/
/////////////////////////////////////////////////////



freeze;

declare verbose AlgEtQIdlWkClasses,3;

///# Weak equivalence classes
/// Let $R$ be an order $I$ and $J$ be fractional $R$-ideals in an étale algebra $A$ over $\mathbb{Q}$. Then the following statements are equivalent:
/// - $I_{\mathfrak{p}} \simeq J_{\mathfrak{p}}$ as $R_{\mathfrak{p}}$- modules for every maximal ideal $\mathfrak{p}$ of $R$.
/// - $I_p \simeq J_p$ as $R_p$- modules for every rational prime $p$. 
/// - $1 \in (I:J)(J:I)$.
/// - $I$ and $J$ have the same multiplicator ring, say $S$, and $(I:J)$ is an invertible fractional $S$-ideal.
/// - $I$ and $J$ have the same multiplicator ring, say $S$, and $(J:I)$ is an invertible fractional $S$-ideal.
/// If the previous conditions hold, one says that $I$ and $J$ are `weakly equivalent`.
/// This is an equivalence relation.
/// 
/// Note that if $I$ is an order, then $J$ is weakly equivalent to $I$ if and only if $J$ is an invertible $I$-ideal. If both $I$ and $J$ are orders then they are weakly equivalent if and only if they are equal.
/// 
/// In the literature, the same notion sometimes goes under the name `local isomorphism`.
/// Reference: Stefano Marseglia, "Local isomorphism classes of fractional ideals of orders in étale algebras", Journal of Algebra 673 (2025), 77-102.

/// Given fractional ideals $I$ and $J$, returns whether they are weakly equivalent. The fractional ideals are no required to be defined over the same order.}
intrinsic IsWeakEquivalent(I::AlgEtQIdl,J::AlgEtQIdl)->BoolElt
{Given fractional ideals I and J, returns whether they are weakly equivalent, that is, if 1 \in (I:J)*(J:I). The fractional ideals are no required to be defined over the same order.}
    S := MultiplicatorRing(I);
    if MultiplicatorRing(J) ne S then
        return false;
    else
        IS:=S!!I;
        JS:=S!!J;
        CIJ:=ColonIdeal(IS,JS);
        CJI:=ColonIdeal(JS,IS);
        //test := OneIdeal(S) eq (CIJ*CJI); //note that this test does not depend on the order of definition of the ideals.
        id:=(CIJ*CJI);
        test:=One(Algebra(I)) in id; //faster!
        return test;
    end if;
end intrinsic;

///ditto
intrinsic IsWeakEquivalent(O1::AlgEtQOrd,O2::AlgEtQOrd)->BoolElt
{Check if the two orders are weakly equivalent, that is equal.}
    return O1 eq O2;
end intrinsic;

///ditto
intrinsic IsWeakEquivalent(O::AlgEtQOrd,J::AlgEtQIdl)->BoolElt
{Checks if the ideal J is weakly equivalent to order O, that is, if J is invertible in O.}
    return IsWeakEquivalent(OneIdeal(O), J);
end intrinsic;

///ditto
intrinsic IsWeakEquivalent(J::AlgEtQIdl,O::AlgEtQOrd)->BoolElt
{Checks if the ideal J is weakly equivalent to order O, that is, if J is invertible in O.}
    return IsWeakEquivalent(OneIdeal(O), J);
end intrinsic;

///ditto
intrinsic IsWeaklyEquivalent(I::AlgEtQIdl,J::AlgEtQIdl)->BoolElt
{Checks if I and J are weakly equivalent, that is, if 1 \in (I:J)*(J:I), or equivalently, if I and J are locally equivalent at all prime of their common multiplicator ring. This function does not require that the ideals are defined over the same order.}
    return IsWeakEquivalent(I,J);
end intrinsic;

///ditto
intrinsic IsWeaklyEquivalent(O1::AlgEtQOrd,O2::AlgEtQOrd)->BoolElt
{Check if the two orders are weakly equivalent, that is equal.}
    return O1 eq O2;
end intrinsic;

///ditto
intrinsic IsWeaklyEquivalent(O::AlgEtQOrd,J::AlgEtQIdl)->BoolElt
{Checks if the ideal J is weakly equivalent to order O, that is, if J is invertible in O.}
    return IsWeakEquivalent(OneIdeal(O), J);
end intrinsic;

///ditto
intrinsic IsWeaklyEquivalent(J::AlgEtQIdl,O::AlgEtQOrd)->BoolElt
{Checks if the ideal J is weakly equivalent to order O, that is, if J is invertible in O.}
    return IsWeakEquivalent(OneIdeal(O), J);
end intrinsic;

/* TESTS

    printf "### Testing WKEq:";
	//AttachSpec("~/packages_github/AlgEt/spec");
	SetAssertions(2);
	_<x>:=PolynomialRing(Integers());
    f:=x^4-100*x^3-100*x^2-100*x-100;
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    assert not IsWeakEquivalent(E,Conductor(E));
    assert not IsWeakEquivalent(OneIdeal(E),Conductor(E));
    assert not IsWeakEquivalent(E,MaximalOrder(K));
    assert IsWeakEquivalent(OneIdeal(MaximalOrder(K)),Conductor(E));
    SetAssertions(1);
    printf " all good!"; 
*/