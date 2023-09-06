/* vim: set syntax=magma :*/

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
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

declare attributes AlgEtQ:TotallyRealSubAlgebra;
declare attributes AlgEtQOrd:TotallyRealUnitGroup;
declare attributes AlgEtQOrd:TotallyRealPositiveUnitGroup;

//////////////////////////
// IsTotallyReal, IsTotallyRealPositive
//////////////////////////

    intrinsic IsTotallyReal(a::AlgEtQElt) -> BoolElt
    {Returns whether a is totally real.}
        return a eq ComplexConjugate(a); 
    end intrinsic;

    intrinsic IsTotallyRealPositive(a::AlgEtQElt) -> BoolElt
    {Returns whether a is totally positive, that is, totally real and with positive image in CC.}
        return IsTotallyReal(a) and forall{ h : h in HomsToC(Parent(a)) | Re(h(a)) gt 0 }; 
    end intrinsic;

//////////////////////////
// Totally Real SubAlgebra
//////////////////////////

    intrinsic TotallyRealSubAlgebra(K::AlgEtQ) -> AlgEtQ,Map
    {Given a CM algebra K, returns the unique totally real subalgebra, with an embedding.}
        if not assigned K`TotallyRealSubAlgebra then
            require HasComplexConjugate(K) : "the algebra does not have CM ";
            a:=PrimitiveElement(K);
            b:=a+ComplexConjugate(a);
            F:=EtaleAlgebra(MinimalPolynomial(b));
            dim:=Dimension(F); 
            bF:=PrimitiveElement(F);
            bF_pows:=[i gt 1 select Self(i-1)*bF else One(F) : i in [1..dim]];
            b_pows:=[i gt 1 select Self(i-1)*b else One(K) : i in [1..dim]];
            assert Dimension(K)/Dimension(F) eq 2;
            FtoK:=map<F->K | x:->SumOfProducts(b_pows,AbsoluteCoordinates([x],bF_pows)[1]) >;
            //FtoK:=hom<F->K | [&+[b^i*coord[j][i+1] : i in [0..dim-1]] : j in [1..dim] ] >;
            assert FtoK(PrimitiveElement(F)) eq b; 
            K`TotallyRealSubAlgebra:=<F,FtoK>;
        end if;
        return K`TotallyRealSubAlgebra[1],K`TotallyRealSubAlgebra[2];
    end intrinsic;

//////////////////////////
// Totally Real Units
//////////////////////////

    intrinsic TotallyRealUnitGroup(S::AlgEtQOrd) -> Grp
    {Given an order S in a CM étale algebra A returns the groups of totally real units of S, as a subgroup of S^*.}
        if not assigned S`TotallyRealUnitGroup then
            K:=Algebra(S);
            F,FtoK:=TotallyRealSubAlgebra(K);
            OF:=MaximalOrder(F); // computing OF is not optimal, and possibly avoidable
                                 // on the other hand, to compu te US, one needs OK
            UF,uF:=UnitGroup(OF); // uF:UF->F
            UK,uK:=UnitGroup(MaximalOrder(K));
            UFtoUK:=hom< UF->UK | [ (FtoK(uF(UF.i))@@uK) : i in [1..Ngens(UF)]]>;
            US,uS:=UnitGroup(S);
            US_in_UK:=sub<UK| [(US.i)@uS@@uK : i in [1..Ngens(US)]]>;

            S_tot_real_in_UK:=Image(UFtoUK) meet US_in_UK;
            S_tot_real_in_US:=sub<US| [ (uK(S_tot_real_in_UK.i))@@uS : i in [1..Ngens(S_tot_real_in_UK)] ]>;
            S`TotallyRealUnitGroup:=S_tot_real_in_US;
        end if;
        return S`TotallyRealUnitGroup;
    end intrinsic;

//////////////////////////
// Totally Real Positive Units
//////////////////////////

    debug_test:=function(S,S_tot_pos_in_US) 
    // this test is quite costly!
    // it loops over the subgroups of S*/<v*bar(v) : v in S^*> meet S^* to find the totally positive units
        K:=Algebra(S);
        zbS:=ZBasis(S);
        T:=Order(zbS cat [ ComplexConjugate(z) : z in zbS ]); // T=S*bar(S)
        UT,uT:=UnitGroup(T); //uT:UT->T
        US, uS := UnitGroup(S); //uS:US->S
        gensUinS:=[ uS(US.i) : i in [1..Ngens(US)]];
        USUSb:=sub< UT | [ (g*ComplexConjugate(g))@@uT : g in gensUinS ]>; // <v*bar(v) : v in S^*> in UT
        USinUT:=sub<UT | [ g@@uT : g in gensUinS ]>;
        Q,q:=quo< USinUT | USinUT meet USUSb >; // q:=USinUT->Q
                                                // Q = S*/<v bar(v) : v in S*> meet S* in UT
        ker_q:=Kernel(q);
        subs:=Subgroups(Q);
        subs_tp:=[]; // subs of tot pos elements
        for H0 in subs do
            H:=H0`subgroup;
            gensH_inK:=[ uT((Q!H.i)@@q) : i in [1..Ngens(H)]];
            if forall{ g : g in gensH_inK | IsTotallyRealPositive(g)} then
                Append(~subs_tp,H);
            end if;
        end for;
        Htp:=&+subs_tp; // in Q 
        Htp:=Htp@@q; // in USinUT
        gens_Htp_inK:=[Htp.i@uT : i in [1..Ngens(Htp)]];
        assert forall{ g : g in gens_Htp_inK | IsTotallyRealPositive(g) };
        Htp:=sub<US| [g@@uS : g in gens_Htp_inK ]>;
        return Htp eq S_tot_pos_in_US;
    end function;

    intrinsic TotallyRealPositiveUnitGroup(S::AlgEtQOrd) -> Grp
    {Given an order S in a CM étale algebra. Returns the groups of totally positive units of S, as a subgroup of S^*.}
        if not assigned S`TotallyRealPositiveUnitGroup then
            K:=Algebra(S);
            F,FtoK:=TotallyRealSubAlgebra(K);
            OF:=MaximalOrder(F); // computing OF is not optimal, and possibly avoidable
                                 // on the other hand, to compute US, one needs OK
            UF,uF:=UnitGroup(OF); // uF:UF->F
            UK,uK:=UnitGroup(MaximalOrder(K));
            UFtoUK:=hom< UF->UK | [ (FtoK(uF(UF.i))@@uK) : i in [1..Ngens(UF)]]>;
            US,uS:=UnitGroup(S);
            US_in_UK:=sub<UK| [(US.i)@uS@@uK : i in [1..Ngens(US)]]>;

            FF2:=FiniteField(2);
            signs:=Matrix([[ Re(h(uF(UF.i))) gt 0 select (FF2!0) else (FF2!1) : h in HomsToC(F) ] : i in [1..Ngens(UF)]]); 
            sol,ker:=Solution( signs , Vector([FF2!0 : i in [1..Dimension(F)]]));
            UF_tot_pos_gens:=   [ 2*UF.i : i in [1..Ngens(UF)] ] cat
                                [ &+[ sol[i]+k[i] eq 1 select UF.i else 0*UF.i : i in [1..Ngens(UF)] ] : k in Basis(ker)] ;
            assert2 debug_test(OF,sub<UF|UF_tot_pos_gens>);
            S_tot_pos_in_UK:=sub<UK| [ g@uF@FtoK@@uK : g in UF_tot_pos_gens ]> meet US_in_UK;
            S_tot_pos_in_US:=sub<US| [ (uK(S_tot_pos_in_UK.i))@@uS : i in [1..Ngens(S_tot_pos_in_UK)] ]>;
            assert2 debug_test(S,S_tot_pos_in_US);
            S`TotallyRealPositiveUnitGroup:=S_tot_pos_in_US;
        end if;
        return S`TotallyRealPositiveUnitGroup;
    end intrinsic;

/* TESTS

    printf "### Testing TotRealPos:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    PP<x>:=PolynomialRing(Integers());
    SetAssertions(2);

    f:=x^8+16;
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    R:=Order([F,2/F]);
    oo:=FindOverOrders(R);
    for iS in [1..#oo] do
        S:=oo[iS];
        _:=TotallyRealUnitGroup(S);
        _:=TotallyRealPositiveUnitGroup(S);
    end for;
    SetAssertions(1);
    printf " all good!\n"; 

*/
