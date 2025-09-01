/////////////////////////////////////////////////////
// Stefano Marseglia, stefano.marseglia89@gmail.com
// https://stmar89.github.io/index.html
// 
// Distributed under the terms of the GNU Lesser General Public License (L-GPL)
//      http://www.gnu.org/licenses/
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 3.0 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301  USA
// 
// Copyright 2023, S. Marseglia
/////////////////////////////////////////////////////

freeze;

declare verbose ComplexConj, 1;

declare attributes AlgEtQ : HasComplexConjugate;
declare attributes AlgEtQOrd : IsConjugateStable;
declare attributes AlgEtQIdl : IsConjugateStable;

///# Complex Conjugation and CM-étale algebras
/// Let $A$ be an étale algebra over $\mathbb{Q}$ with components $K_1\times\cdots\times K_n$.
/// We say that $A$ is a `CM-étale algebra` if every component $K_i$ is a CM-field, that is, $K_i$ has involution that acts as applying complex conjugation after applying any homomorphism to the complex numbers. If $A$ is a CM-étale algebra, then it has an involution with the same property. For this reason, we call this involution `complex conjugation`.
/// Given an element of $A$, an order or a fractional ideal in $A$, we say that it is `conjugate stable` if it equals it complex conjugate.

intrinsic HasComplexConjugate(A::AlgEtQ) -> BoolElt
{Returns if the algebra is the product of CM fields.}
    if not assigned A`HasComplexConjugate then
        A`HasComplexConjugate:=forall{L : L in A`Components[1] | HasComplexConjugate(L)};
    end if;
	return A`HasComplexConjugate;
end intrinsic;

intrinsic ComplexConjugate(x::AlgEtQElt) -> AlgEtQElt
{Returns the complex conjugate of the argument.}
	A:=Parent(x);
	require HasComplexConjugate(A) : "it is not a product of CM fields";
    return A!<ComplexConjugate(xi) : xi in Components(x)>;
end intrinsic;

/// Given an order $O$ in a CM-étale algebra, returns wheter $O$ is conjugate stable and the complex conjugate.
intrinsic IsConjugateStable(O::AlgEtQOrd) -> BoolElt,AlgEtQOrd
{Given an order O in a CM-étale algebra, returns wheter O is conjugate stable and the complex conjugate.}
    if not assigned O`IsConjugateStable then
        A:=Algebra(O);
        Ob:=Order([ ComplexConjugate(x) : x in ZBasis(O) ] : CheckIsKnownOrder:=false);
        is_stable:=Ob eq O;
        if is_stable then
            O`IsConjugateStable:=<true,O>;
        else
            IsKnownOrder(~Ob);
            Ob`IsConjugateStable:=<false,O>;
            O`IsConjugateStable:=<false,Ob>;
        end if;
    end if;
	return Explode(O`IsConjugateStable);
end intrinsic;

/// Given an order $O$ in a CM-étale algebra, returns the complex conjugate of $O$.
intrinsic ComplexConjugate(O::AlgEtQOrd) -> AlgEtQOrd
{Given an order O in a CM-étale algebra, returns the complex conjugate of O.}
	A:=Algebra(O);
    _,Ob:=IsConjugateStable(O); //if stable Ob = O, to preserve atttributes!
	return Ob;
end intrinsic;

/// Given a fractional ideal $I$ in a CM-étale algebra, returns wheter $I$ is conjugate stable and the complex conjugate. If the order of $I$ is not conjugate stable, then the second output will be defined over the complex conjugate of the order.
intrinsic IsConjugateStable(I::AlgEtQIdl) -> BoolElt,AlgEtQIdl
{Given a fractional ideal I in a CM-étale algebra, returns wheter I is conjugate stable and the complex conjugate. If the order of I is not conjugate stable, then the second output will be defined over the complex conjugate of the order.}
    if not assigned I`IsConjugateStable then
        O:=Order(I);
        A:=Algebra(O);
        is_stableO,Ob:=IsConjugateStable(O);
        Ib:=Ideal(Ob,[ ComplexConjugate(x) : x in ZBasis(I)]);
        if not is_stableO then
            is_stable:=false; //can't compare I and Ib sicne they are defined over different orders
        else
            is_stable:=I eq Ib;
        end if;
        if is_stable then
            I`IsConjugateStable:=<true,I>;
        else
            I`IsConjugateStable:=<false,Ib>;
        end if;
    end if;
	return Explode(I`IsConjugateStable);
end intrinsic;

/// Given a fractional ideal $I$ in a CM-étale algebra, returns the complex conjugate of $I$. If the order of $I$ is not conjugate stable, then the second output will be defined over the complex conjugate of the order.
intrinsic ComplexConjugate(I::AlgEtQIdl) -> AlgEtQIdl
{Given a fractional ideal I in a CM-étale algebra, returns the complex conjugate of I. If the order of I is not conjugate stable, then the second output will be defined over the complex conjugate of the order.}
    _,Ib:=IsConjugateStable(I); //if stable Ib = I, to preserve attributes!
	return Ib;
end intrinsic;

/* TESTS

    printf "### Testing Complex Conjugation:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    _<x>:=PolynomialRing(Integers());
    f:=x^4 + 6*x^2 + 25;
    K:=EtaleAlgebra(f);
    pi:=PrimitiveElement(K);
    pib:=ComplexConjugate(pi);
    assert pib eq 5/pi;
    R:=Order([pi,pib]);
    oo:=FindOverOrders(R);
    O:=MaximalOrder(K);
    S:=[ S : S in oo | Index(O,S) eq 8 ][1]; //there is only one order of index 8, so conjugate stable
    assert IsConjugateStable(R);
    assert IsConjugateStable(TraceDualIdeal(R));
    assert IsConjugateStable(O);
    assert IsConjugateStable(TraceDualIdeal(O));
    assert IsConjugateStable(S);
    assert IsConjugateStable(TraceDualIdeal(S));
    assert not IsConjugateStable(EquationOrder(K));
    printf ".";
    printf " all good!\n";

*/
