/* vim: set syntax=magma :*/

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

declare attributes AlgEtQOrd:   ICM,
                                ICM_bar;

intrinsic ICM_bar(S::AlgEtQOrd : GRH:=false ) -> SeqEnum
{Returns the ideal classes of fractional S-ideals having MultiplicatorRing equal to S. This is the same as the orbit of the action of PicardGroup(S) on WKICM_bar(S).}
    if not assigned S`ICM_bar then
        seqWKS_bar:=WKICM_bar(S);
        GS,gS:=PicardGroup(S : GRH:=GRH );
        repS:=[gS(x) : x in GS];
        ICM_barS := &cat[[I*J : I in seqWKS_bar] : J in repS];
        for iI in [1..#ICM_barS] do
            I:=ICM_barS[iI];
            I`MultiplicatorRing:=S;
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
        _:=WKICM(S); //to populate all attributes T`WKICM_bar
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




