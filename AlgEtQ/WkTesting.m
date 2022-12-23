/* vim: set syntax=magma :*/

//freeze;

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
{Checks if I and J are weakly equivalent 1 \in (I:J)*(J:I). This function does not require that the ideals are defined over the same order.}
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
{Checks if the second argument is weakly equivalent to the first argument.}
    return IsWeakEquivalent(OneIdeal(O), J);
end intrinsic;

intrinsic IsWeakEquivalent(J::AlgEtQIdl,O::AlgEtQOrd)->BoolElt
{Checks if the second argument is weakly equivalent to the first argument.}
    return IsWeakEquivalent(OneIdeal(O), J);
end intrinsic;

intrinsic IsGorenstein(O::AlgEtQOrd)->BoolElt
{Checks if the order O is Gorenstein.}
    if not assigned O`IsGorenstein then
        if assigned O`IsMaximal and O`IsMaximal then
            O`IsGorenstein:=true;
        else
            T:=TraceDualIdeal(O);
            O`IsGorenstein:=IsInvertible(T);
        end if;
    end if;
    return O`IsGorenstein;
end intrinsic



/* TEST

    printf "### Testing WKICM:";
	AttachSpec("~/packages_github/AlgEt/spec");
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
