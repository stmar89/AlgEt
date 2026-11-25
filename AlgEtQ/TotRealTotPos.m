/////////////////////////////////////////////////////
// Copyright 2025.
// Stefano Marseglia, stefano.marseglia89@gmail.com
// https://stmar89.github.io/index.html
//
// Distributed under the terms of the CC-BY 4.0 licence.
// https://creativecommons.org/licenses/by/4.0/
/////////////////////////////////////////////////////



freeze;

declare attributes AlgEtQ:TotallyRealSubAlgebra;
declare attributes AlgEtQOrd:TotallyRealUnitGroup;
declare attributes AlgEtQOrd:TotallyRealPositiveUnitGroup;

import "PicardGroup.m" : conductor_sequence;

units_of_S_in_units_of_O:=function(S)
// Given an order S, returns US,uS:=UnitGroup(S) and the unit group UK,uK of the maximal order, such
// that US is a subgroup of UK. No unit of S is computed: for a non-maximal S, UnitGroup(S) is the
// kernel of a map from UK (see conductor_sequence in PicardGroup.m). This matters when the
// index of US in UK is huge: then the units of S have huge exponents on the generators of UK, and
// they cannot be written down as elements of the algebra.
    US,uS:=UnitGroup(S);
    if IsMaximal(S) then
        return US,uS,US,uS;
    end if;
    data:=conductor_sequence(S,false); // already computed by UnitGroup(S)
    assert IsIdentical(US,data`K); // US is the subgroup data`K of data`UO
    return US,uS,data`UO,data`uO;
end function;

totally_real_units_data:=function(S)
// Given an order S in a CM algebra K, returns the totally real subalgebra F of K with its embedding
// FtoK, the unit group UF,uF of the maximal order of F, US:=UnitGroup(S), and the unit group UK,uK
// of the maximal order of K, of which US is a subgroup (see units_of_S_in_units_of_O). Shared by
// TotallyRealUnitGroup and TotallyRealPositiveUnitGroup.
    K:=Algebra(S);
    F,FtoK:=TotallyRealSubAlgebra(K);
    OF:=MaximalOrder(F); // computing OF is not optimal, and possibly avoidable
                         // on the other hand, to compute US, one needs OK
    UF,uF:=UnitGroup(OF); // uF:UF->F
    US,_,UK,uK:=units_of_S_in_units_of_O(S); // US is a subgroup of UK
    return F,FtoK,UF,uF,US,UK,uK;
end function;

meet_in_US:=function(H, US, UK)
// Given subgroups H and US of UK, returns H meet US as a subgroup of US. The intersection is taken
// in UK, and written in US by coercion: no unit of S is computed.
    H_meet_US:=H meet US;
    return sub<US| [ US!(UK!H_meet_US.i) : i in [1..Ngens(H_meet_US)] ]>;
end function;

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
        FtoK:=map<F->K | x:->DotProduct(b_pows,AbsoluteCoordinates([x],bF_pows)[1]) >;
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
{Given an order S in a CM étale algebra A returns the group of totally real units of S, as a subgroup of S^*.}
    if not assigned S`TotallyRealUnitGroup then
        F,FtoK,UF,uF,US,UK,uK:=totally_real_units_data(S);
        UFtoUK:=hom< UF->UK | [ (FtoK(uF(UF.i))@@uK) : i in [1..Ngens(UF)]]>;
        S`TotallyRealUnitGroup:=meet_in_US(Image(UFtoUK),US,UK);
    end if;
    return S`TotallyRealUnitGroup;
end intrinsic;

//////////////////////////
// Totally Real Positive Units
//////////////////////////

debug_test:=function(S,S_tot_pos_in_US)
// this test is quite costly!
// it loops over the subgroups of S*/<v*bar(v) : v in S^*> meet S^* to find the totally positive
// units
    zbS:=ZBasis(S);
    T:=Order(zbS cat [ ComplexConjugate(z) : z in zbS ]); // T=S*bar(S)
    UT,uT:=UnitGroup(T); //uT:UT->T
    US, uS := UnitGroup(S); //uS:US->S
    gensUinS:=[ uS(US.i) : i in [1..Ngens(US)]];
    // <v*bar(v) : v in S^*> in UT
    USUSb:=sub< UT | [ (g*ComplexConjugate(g))@@uT : g in gensUinS ]>;
    USinUT:=sub<UT | [ g@@uT : g in gensUinS ]>;
    Q,q:=quo< USinUT | USinUT meet USUSb >; // q:=USinUT->Q
                                            // Q = S*/<v bar(v) : v in S*> meet S* in UT
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
{Given an order S in a CM étale algebra. Returns the group of totally positive totally real units of S, as a subgroup of S^*.}
    if not assigned S`TotallyRealPositiveUnitGroup then
        F,FtoK,UF,uF,US,UK,uK:=totally_real_units_data(S);
        FF2:=FiniteField(2);
        signs:=Matrix([[ Re(h(uF(UF.i))) gt 0 select (FF2!0) else (FF2!1) : h in HomsToC(F) ] : i in [1..Ngens(UF)]]);
        sol,ker:=Solution( signs , Vector([FF2!0 : i in [1..Dimension(F)]]));
        UF_tot_pos_gens:=   [ 2*UF.i : i in [1..Ngens(UF)] ] cat
                            [ &+[ sol[i]+k[i] eq 1 select UF.i else 0*UF.i : i in [1..Ngens(UF)] ] : k in Basis(ker)] ;
        assert2 debug_test(MaximalOrder(F),sub<UF|UF_tot_pos_gens>);
        S_tot_pos_in_US:=meet_in_US(sub<UK| [ g@uF@FtoK@@uK : g in UF_tot_pos_gens ]>,US,UK);
        assert2 debug_test(S,S_tot_pos_in_US);
        S`TotallyRealPositiveUnitGroup:=S_tot_pos_in_US;
    end if;
    return S`TotallyRealPositiveUnitGroup;
end intrinsic;

/* TESTS

    printf "### Testing TotRealPos:";
    //AttachSpec("../spec");
    PP<x>:=PolynomialRing(Integers());
    SetAssertions(2);

    // The method used before 2026-11: compute the units of S and take their discrete logarithm in
    // the unit group of the maximal order. Now US is used directly as a subgroup of it. The two
    // must give the same subgroups.
    old_tot_real:=function(S)
        K:=Algebra(S);
        F,FtoK:=TotallyRealSubAlgebra(K);
        UF,uF:=UnitGroup(MaximalOrder(F));
        UK,uK:=UnitGroup(MaximalOrder(K));
        UFtoUK:=hom< UF->UK | [ (FtoK(uF(UF.i))@@uK) : i in [1..Ngens(UF)]]>;
        US,uS:=UnitGroup(S);
        US_in_UK:=sub<UK| [(US.i)@uS@@uK : i in [1..Ngens(US)]]>;
        T:=Image(UFtoUK) meet US_in_UK;
        return sub<US| [ (uK(T.i))@@uS : i in [1..Ngens(T)] ]>;
    end function;
    old_tot_pos:=function(S)
        K:=Algebra(S);
        F,FtoK:=TotallyRealSubAlgebra(K);
        UF,uF:=UnitGroup(MaximalOrder(F));
        UK,uK:=UnitGroup(MaximalOrder(K));
        US,uS:=UnitGroup(S);
        US_in_UK:=sub<UK| [(US.i)@uS@@uK : i in [1..Ngens(US)]]>;
        FF2:=FiniteField(2);
        signs:=Matrix([[ Re(h(uF(UF.i))) gt 0 select (FF2!0) else (FF2!1) : h in HomsToC(F) ] : i in [1..Ngens(UF)]]);
        sol,ker:=Solution( signs , Vector([FF2!0 : i in [1..Dimension(F)]]));
        gens:=[ 2*UF.i : i in [1..Ngens(UF)] ] cat
              [ &+[ sol[i]+k[i] eq 1 select UF.i else 0*UF.i : i in [1..Ngens(UF)] ] : k in Basis(ker)];
        T:=sub<UK| [ g@uF@FtoK@@uK : g in gens ]> meet US_in_UK;
        return sub<US| [ (uK(T.i))@@uS : i in [1..Ngens(T)] ]>;
    end function;

    f:=x^8+16;
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    R:=Order([F,2/F]);
    oo:=FindOverOrders(R);
    for iS in [1..#oo] do
        S:=oo[iS];
        tr:=TotallyRealUnitGroup(S);
        tp:=TotallyRealPositiveUnitGroup(S);
        assert tr eq old_tot_real(S);
        assert tp eq old_tot_pos(S);
        // the generators are totally real, and totally positive
        US,uS:=UnitGroup(S);
        assert forall{ g : g in Generators(tr) | IsTotallyReal(uS(US!g)) };
        assert forall{ g : g in Generators(tp) | IsTotallyRealPositive(uS(US!g)) };
    end for;

    // An order whose unit group has huge index in the one of the maximal order: [O^*:E^*] has 20
    // digits. Its units of infinite order cannot be written down as elements, so the old method
    // fails, while the subgroups are computed inside the unit group of the maximal order. With
    // assertions at level 1, since debug_test (assert2) computes the units of E.
    SetAssertions(1);
    A:=EtaleAlgebra([NumberField(x^4+x^2+529),NumberField(x^8+16)]);
    E:=EquationOrder(A);
    assert Invariants(TotallyRealUnitGroup(E)) eq [2,0,0,0,0];
    assert Invariants(TotallyRealPositiveUnitGroup(E)) eq [0,0,0,0];
    SetAssertions(1);
    printf " all good!";
*/
