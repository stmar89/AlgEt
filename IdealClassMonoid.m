/* vim: set syntax=magma :*/

//freeze;

/////////////////////////////////////////////////////
// Ideal class monoid and weak equivalence classes for orders in Etale Q algebras
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare attributes AlgEtOrd:ICM,
                            ICM_bar;

intrinsic ICM_bar(S::AlgEtOrd : GRH:=false ) -> SeqEnum
{returns the ideal classes of the order S having S as MultiplicatorRing, that is the orbits of the action of PicardGroup(S) on WKICM_bar(S)}
    if not assigned S`ICM_bar then
        seqWKS_bar:=WKICM_bar(S);
        GS,gS:=PicardGroup(S : GRH:=GRH );
        repS:=[gS(x) : x in GS];
        ICM_barS := &cat[[I*J : I in seqWKS_bar] : J in repS];
        assert2 forall{J : J in ICM_barS | MultiplicatorRing(J) eq S};
        assert forall{J : J in ICM_barS | Order(J) eq S};
        S`ICM_bar:=ICM_barS;
    end if;
    return S`ICM_bar;
end intrinsic;

intrinsic ICM(S::AlgEtOrd : GRH:=false ) -> SeqEnum
{returns the ideal class monoid of the order, that is a set of representatives for the isomorphism classes of the fractiona ideals}
    if not assigned S`ICM then
        seqOO:=FindOverOrders(S);
        seqICM:=[];
        for T in seqOO do
            ICM_barT := [(S!!I) : I in ICM_bar(T: GRH:=GRH )];
            seqICM:=seqICM cat ICM_barT;
        end for;
        assert forall{I: I in seqICM | Order(I) eq S};
        S`ICM:=seqICM;
    end if;
    return S`ICM;
end intrinsic;

/*TEST


*/




