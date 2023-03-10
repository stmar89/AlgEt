/* vim: set syntax=magma :*/

freeze;

/////////////////////////////////////////////////////
// Ideal class monoid and weak equivalence classes for orders in Etale Q algebras
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare attributes AlgEtQOrd:   ICM,
                                ICM_bar;

intrinsic ICM_bar(S::AlgEtQOrd : GRH:=false ) -> SeqEnum
{Returns the ideal classes of fractional S-ideals having MultiplicatorRing equal to S. This is the same as the orbit of the action of PicardGroup(S) on WKICM_bar(S).}
    if not assigned S`ICM_bar then
        seqWKS_bar:=WKICM_bar(S);
        GS,gS:=PicardGroup(S : GRH:=GRH );
        repS:=[gS(x) : x in GS];
        ICM_barS := &cat[[I*J : I in seqWKS_bar] : J in repS];
        for I in ICM_barS do
            ZBasisLLL(I);
        end for;
        assert2 forall{J : J in ICM_barS | MultiplicatorRing(J) eq S};
        assert forall{J : J in ICM_barS | Order(J) eq S};
        S`ICM_bar:=ICM_barS;
    end if;
    return S`ICM_bar;
end intrinsic;

intrinsic ICM(S::AlgEtQOrd : GRH:=false ) -> SeqEnum
{Returns the ideal class monoid of the order, that is, a set of representatives for the isomorphism classes of the fractional S-ideals.}
    if not assigned S`ICM then
        _:=WKICM(S); // This is to populate the attribute WKICM_bar. Currently it is superfluous.
                     // It will make the computation faster once the recursive method for WKICM is included.
        seqOO:=OverOrders(S);
        seqICM:=[];
        for T in seqOO do
            ICM_barT := [(S!!I) : I in ICM_bar(T: GRH:=GRH )];
            seqICM:=seqICM cat ICM_barT;
        end for;
        assert forall{I: I in seqICM | Order(I) eq S};
        for I in seqICM do
            ZBasisLLL(I);
        end for;
        S`ICM:=seqICM;
    end if;
    return S`ICM;
end intrinsic;

/* TESTS

    printf "### Testing ICM:";
	SetAssertions(2);
    SetClassGroupBounds("GRH");
	_<x>:=PolynomialRing(Integers());
    polys:=[
        x^3-100*x^2-100*x-100,
        x^4+291*x^3-988*x^2-1000*x-1000,
        x^3+31*x^2+43*x+77
        ];
    for f in polys do
        K:=EtaleAlgebra(f);
        E:=EquationOrder(K);
        _:=ICM(E);
    end for;
    SetAssertions(1);
    printf " all good!\n";

*/




