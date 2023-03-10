/* vim: set syntax=magma :*/

freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose AlgEtQIdlWkClasses,3;

/*TODO:

*/

//----------
// Weak Equiv testing
//----------

intrinsic IsWeakEquivalent(I::AlgEtQIdl,J::AlgEtQIdl)->BoolElt
{Checks if I and J are weakly equivalent, that is, if 1 \in (I:J)*(J:I), or equivalently, if I and J are locally equivalent at all prime of their common multiplicator ring. This function does not require that the ideals are defined over the same order.}
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

intrinsic IsWeakEquivalent(O1::AlgEtQOrd,O2::AlgEtQOrd)->BoolElt
{Check if the two orders are weakly equivalent, that is equal.}
    return O1 eq O2;
end intrinsic;

intrinsic IsWeakEquivalent(O::AlgEtQOrd,J::AlgEtQIdl)->BoolElt
{Checks if the ideal J is weakly equivalent to order O, that is, if J is invertible in O.}
    return IsWeakEquivalent(OneIdeal(O), J);
end intrinsic;

intrinsic IsWeakEquivalent(J::AlgEtQIdl,O::AlgEtQOrd)->BoolElt
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
    printf " all good!\n"; 
    

*/
