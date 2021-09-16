/* vim: set syntax=magma :*/

//freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose AlgEtIdlWkClasses,3;

/*TODO:

*/

//----------
// Weak Equiv testing
//----------

intrinsic IsWeakEquivalent(I::AlgEtIdl,J::AlgEtIdl)->BoolElt
{ checks if 1 \in (I:J)*(J:I). This function does not require that the ideals are defined over the same order. }
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

intrinsic IsWeakEquivalent(O1::AlgEtOrd,O2::AlgEtOrd)->BoolElt
{ check if the two orders are weakly equivalent, that is equal }
    return O1 eq O2;
end intrinsic;

intrinsic IsWeakEquivalent(O::AlgEtOrd,J::AlgEtIdl)->BoolElt
{ checks if the second argument is weakly equivalent to the first argument }
    return IsWeakEquivalent(OneIdeal(O), J);
end intrinsic;

intrinsic IsWeakEquivalent(J::AlgEtIdl,O::AlgEtOrd)->BoolElt
{ checks if the second argument is weakly equivalent to the first argument }
    return IsWeakEquivalent(OneIdeal(O), J);
end intrinsic;

intrinsic IsInvertible(I::AlgEtIdl) ->BoolElt,AlgEtIdl
{ checks if the ideal I is invertible in its order of definition O, and the colon ideal (O:I) }
    if not assigned I`IsInvertible then
        O:=Order(I);
        COI:=ColonIdeal(O,I);
        I`IsInvertible:=<O eq I*COI, COI>;
    end if;
    return Explode(I`IsInvertible);
end intrinsic;



/* TEST

*/
