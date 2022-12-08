/* vim: set syntax=magma :*/

freeze;

/////////////////////////////////////////////////////
// Complex Conjugation for AlgEt
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose ComplexConj, 1;

declare attributes AlgEt : HasComplexConjugate;
declare attributes AlgEtOrd : IsConjugateStable;
declare attributes AlgEtIdl : IsConjugateStable;

intrinsic HasComplexConjugate(A::AlgEt) -> BoolElt
{Returns if the algebra is the product of CM fields}
    if not assigned A`HasComplexConjugate then
        A`HasComplexConjugate:=forall{L : L in A`Components[1] | HasComplexConjugate(L)};
    end if;
	return A`HasComplexConjugate;
end intrinsic;

intrinsic ComplexConjugate(x::AlgEtElt) -> AlgEtElt
{If A is a product of CM fields, it returns the complex conjugate of the argument}
	A:=Parent(x);
	require HasComplexConjugate(A) : "it is not a product of CM fields";
    return A!<ComplexConjugate(xi) : xi in Components(x)>;
end intrinsic;

intrinsic IsConjugateStable(O::AlgEtOrd) -> BoolElt,AlgEtOrd
{Returns wheter O is conjugate stable and the complex conjugate.}
    if not assigned O`IsConjugateStable then
        A:=Algebra(O);
        Ob:=Order([ ComplexConjugate(x) : x in ZBasis(O) ]);
        is_stable:=Ob eq O;
        if is_stable then
            O`IsConjugateStable:=<true,O>;
        else
            O`IsConjugateStable:=<false,Ob>;
        end if;
    end if;
	return Explode(O`IsConjugateStable);
end intrinsic;

intrinsic ComplexConjugate(O::AlgEtOrd) -> AlgEtOrd
{It returns the complex conjugate of the argument.}
	A:=Algebra(O);
    _,Ob:=IsConjugateStable(O); //if stable Ob = O, to preserve atttributes!
	return Ob;
end intrinsic;

intrinsic IsConjugateStable(I::AlgEtIdl) -> BoolElt,AlgEtIdl
{Returns wheter O is conjugate stable and the complex conjugate.}
    if not assigned I`IsConjugateStable then
        O:=Order(I);
        A:=Algebra(O);
        is_stableO,Ob:=IsConjugateStable(O);
        Ib:=Ideal(Ob,[ ComplexConjugate(x) : x in ZBasis(I)]);
        if not is_stableO then
            is_stableI:=false; //can't compare I and Ib sicne they are defined over different orders
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

intrinsic ComplexConjugate(I::AlgEtIdl) -> AlgEtIdl
{if A is a product of CM fields, it returns the complex conjugate of the argument}
    _,Ib:=IsConjugateStable(I); //if stable Ib = I, to preserve atttributes!
	return Ib;
end intrinsic;

/*
//TESTS

    AttachSpec("~/packages_github/AlgEt/spec");
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


    wkS:=WKICM_bar(S);
    #wkS;
    S eq Order(ZBasis(2*O));

    for I in wkS do
        GI:=Quotient(I,(1-pi)*I);
        Iv:=TraceDualIdeal(ComplexConjugate(I));
        GIv:=Quotient(Iv,(1-pi)*Iv);
        IsIsomorphic(I,Iv);
        ElementaryDivisors(GI),
        ElementaryDivisors(GIv);
        "\n";
    end for;

*/
