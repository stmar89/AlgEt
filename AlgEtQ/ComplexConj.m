/* vim: set syntax=magma :*/

freeze;

/////////////////////////////////////////////////////
// Complex Conjugation for AlgEtQ
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose ComplexConj, 1;

declare attributes AlgEtQ : HasComplexConjugate;
declare attributes AlgEtQOrd : IsConjugateStable;
declare attributes AlgEtQIdl : IsConjugateStable;

intrinsic HasComplexConjugate(A::AlgEtQ) -> BoolElt
{Returns if the algebra is the product of CM fields.}
    if not assigned A`HasComplexConjugate then
        A`HasComplexConjugate:=forall{L : L in A`Components[1] | HasComplexConjugate(L)};
    end if;
	return A`HasComplexConjugate;
end intrinsic;

intrinsic ComplexConjugate(x::AlgEtQElt) -> AlgEtQElt
{If A is a product of CM fields, it returns the complex conjugate of the argument.}
	A:=Parent(x);
	require HasComplexConjugate(A) : "it is not a product of CM fields";
    return A!<ComplexConjugate(xi) : xi in Components(x)>;
end intrinsic;

intrinsic IsConjugateStable(O::AlgEtQOrd) -> BoolElt,AlgEtQOrd
{Given an order O in a CM-étale algebra, it returns wheter O is conjugate stable and the complex conjugate.}
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

intrinsic ComplexConjugate(O::AlgEtQOrd) -> AlgEtQOrd
{Given an order O in a CM-étale algebra, it returns the complex conjugate of the argument.}
	A:=Algebra(O);
    _,Ob:=IsConjugateStable(O); //if stable Ob = O, to preserve atttributes!
	return Ob;
end intrinsic;

intrinsic IsConjugateStable(I::AlgEtQIdl) -> BoolElt,AlgEtQIdl
{Given a fractional ideal I in a CM-étale algebra, it returns wheter I is conjugate stable and the complex conjugate. Note: if the order of I is not conjugate stable, then the second output will be defined over the complex conjugate of the order.}
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

intrinsic ComplexConjugate(I::AlgEtQIdl) -> AlgEtQIdl
{If A is a product of CM fields, it returns the complex conjugate of the fractional ideal I. Note: if the order of I is not conjugate stable, then the output will be defined over the complex conjugate of the order.}
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
