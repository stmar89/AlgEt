/////////////////////////////////////////////////////
// Copyright 2025.
// Stefano Marseglia, stefano.marseglia89@gmail.com
// https://stmar89.github.io/index.html
//
// Distributed under the terms of the CC-BY 4.0 licence.
// https://creativecommons.org/licenses/by/4.0/
/////////////////////////////////////////////////////



freeze;

declare verbose AlgEtQTraceNorm, 3;

///# Trace and norm
/// Let $A$ be an étale algebra over $\mathbb{Q}$, with components $K_1\times\cdots\times K_n$.
/// We define the `(absolute) trace` on $A$ as the additive map $\mathrm{Tr_{A/\mathbb{Q}}}\colon A\to \mathbb{Q}$ that sends an element $a\in A$ to $\sum_{i=1}^n \mathrm{Tr}_{K_i/\mathbb{Q}}(a)$.
/// Let $m_a$ be the matrix representing the multipliction-by-$a$ on $A$ with respect to any basis of $A$ over $\mathbb{Q}$. Then $\mathrm{Tr}_{A/\mathbb{Q}}(a)$ equals the trace of $m_a$.
///
/// We define the `(absolute) norm` on $A$ as the multiplicative map $\mathrm{N}_{A/\mathbb{Q}}\colon A \to \mathbb{Q}$ by sending a unit $a \in A$ to $\prod_{i=1}^n \mathrm{N}_{K_i/\mathbb{Q}}(a)$ and every zero-divisor to $0$.
/// We have $N_{A/\mathbb{Q}}(a)$ equals the determinant of the matrix $m_a$.

//------------
// Trace and Norm
//------------

intrinsic Trace(x::AlgEtQElt) -> FldRatElt
{Returns the (absolute) trace of the element.}
    return AbsoluteTrace(x);
end intrinsic;

///ditto
intrinsic AbsoluteTrace(x::AlgEtQElt) -> FldRatElt
{ " } // "
    return &+[AbsoluteTrace(y) : y in Components(x)];
end intrinsic;

intrinsic Norm(x::AlgEtQElt) -> FldRatElt
{Returns the (absolute) norm of the element.}
    return AbsoluteNorm(x);
end intrinsic;

///ditto
intrinsic AbsoluteNorm(x::AlgEtQElt) -> FldRatElt
{ " } // "
    return &*[AbsoluteNorm(y) : y in Components(x)];
end intrinsic;


///## Trace dual ideals
/// Let $I$ be an order or a fractional ideal in an étale algebra $A$ over $\mathbb{Q}$.
/// We defined the `trace dual ideal` of $I$ as $I^t=\{ a\in A : \mathrm{Tr}_{A/\mathbb{Q}}(a\cdot I) \subseteq \mathbb{Z} \}$.
/// For fractional ideals $I$ and $J$ and a unit $a\in A$, we have:
/// - if $I \subseteq J$ then $J^t \subseteq I^t$ and $\#(J/I) = \#(I^t/J^t)$;
/// - $(aI)^t = \frac{1}{a}I^t$;
/// - $(I+J)^t = I^t \cap J^t$;
/// - $(I\cap J)^t = I^t + J^t$;
/// - $(I:J)^t = I^t\cdot J$.

//------------
// Trace dual ideal
//------------

intrinsic TraceDualIdeal(I::AlgEtQIdl) -> AlgEtQIdl
{Returns the trace dual ideal of the given fractional ideal.}
    if not assigned I`TraceDualIdeal then
        A:=Algebra(I);
        require PrimeField(A) eq BaseField(A) : "implementend only for algebras over the prime field";
        S:=Order(I);
        B:=ZBasis(I);
        Nnf:=#Components(A);
        n:=#B;
        Q:=MatrixRing(RationalField(), n)![AbsoluteTrace(B[i]*B[j]): i, j in [1..n] ];
        QQ:=Q^-1;
        //BB:=[A ! (&+[ (QQ[i,j]*B[j]): j in [1..n]]) : i in [1..n]] ; //too many coercions
        B_comp:=[Components(b) : b in B];
        BB:=< [(&+[ (QQ[i,j]*B_comp[j][k]): j in [1..n]]) : i in [1..n]] : k in [1..Nnf]>;
        BB:=[ A ! < BB[k][i] : k in [1..Nnf] > : i in [1..n] ];
        assert2 forall{ i : i,j in [1..n] | AbsoluteTrace( B[i]*BB[j] ) eq KroneckerDelta(i,j) };
        It:=Ideal(S,BB);
        It`ZBasis:=BB; //we know that BB is a ZBasis
        if assigned I`MultiplicatorRing then
        // the multiplicator ring of I is the same of its trace dual.
            It`MultiplicatorRing:=I`MultiplicatorRing;
        end if;
        ZBasisLLL(It);
        I`TraceDualIdeal:=It;
    end if;
    return I`TraceDualIdeal;
end intrinsic;

intrinsic TraceDualIdeal(O::AlgEtQOrd) -> AlgEtQIdl
{Returns the trace dual ideal of the given order.}
    if not assigned O`TraceDualIdeal then
        Ot:=TraceDualIdeal(OneIdeal(O));
        O`TraceDualIdeal := Ot;
    end if;
    return O`TraceDualIdeal;
end intrinsic;

/* TESTS

    printf "### Testing Trace and Norm:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    SetVerbose("AlgEtQTraceNorm",1);

    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    for i in [1..100] do
        a:=Random(A);
        b:=Random(A);
        assert Trace(a)+Trace(b) eq Trace(a+b);
        assert Norm(a)*Norm(b) eq Norm(a*b);
    end for;
    printf " all good!";
*/