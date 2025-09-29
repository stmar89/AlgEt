/////////////////////////////////////////////////////
// Copyright 2025.
// Stefano Marseglia, stefano.marseglia89@gmail.com
// https://stmar89.github.io/index.html
//
// Distributed under the terms of the CC-BY 4.0 licence.
// https://creativecommons.org/licenses/by/4.0/
/////////////////////////////////////////////////////



freeze;

declare verbose AlgEtQPicardGroup, 3;
declare verbose ResidueRingUnits, 3;

declare attributes AlgEtQOrd: PicardGroup,
                              UnitGroup;

declare attributes AlgEtQIdl: residue_class_ring_unit_subgroup_generators,
                              ResidueRingUnits;

///# Picard group and unit group
/// Let $R$ be an order in an étale algebra $A$ over $\mathbb{Q}$ with maximal order $\mathcal{O}_A$.
/// We say that two fractional $R$-ideals $I$ and $J$ are `isomorphic` if there exists a unit $a\in A$ such that $I=aJ$. Observe that this happens if and only if $I$ and $J$ are isomorphic as $R$-modules.
/// We refer to the isomorphism class $[I]$ of $I$ as its `ideal class`.
///
/// The set of ideal classes of invertible fractional $R$-ideals forms a group under the operation induced by ideal multiplication. This group is called the `Picard group` of $R$ and denoted $\mathrm{Pic}(R)$.
///
/// If $K_1\times \ldots \times K_n$ are the components of $A$ then $\mathrm{Pic}(\mathcal{O}_A) = \prod_{i=1}^n\mathrm{Cl}(\mathcal{O}_{K_i})$, where $\mathrm{Cl}(\mathcal{O}_{K_i})$ is the class group of the number field $K_i$.
/// Also, the unit group $\mathcal{O}_A^\times$ of $\mathcal{O}_A$ satisfies $\mathcal{O}_A^\times = \mathcal{O}_{K_1}^\times \times \ldots \times \mathcal{O}_{K_n}^\times$.
///
/// The Picard group and the unit group of the order $R$ can be computed using the well-known exact sequence:
/// ```math
/// 1 \to R^\times \to \mathcal{O}_A^\times \to \dfrac{\left( \mathcal{O}/\mathfrak{f} \right)^\times}{\left( \mathcal{O}/\mathfrak{f} \right)^\times} \to \mathrm{Pic}(R)\to \mathrm{Pic}(\mathcal{O}_A) \to 1,
/// ```
/// where the first, second and third map are the natural maps, the fourth is induced by $a \mapsto (a\mathcal{O}_A \cap R)$ and the last one is induced by the extension map $I \mapsto I\mathcal{O}_A$.
///
/// The following intrinsics allow to compute all the elements of the sequence.
/// If one wants to perform the computation under the GRH bound one can either use the appropriate parameters of use `SetClassGroupBounds("GRH")` as for number fields.


/// Given an order $S$ and a fractional $S$-ideal $I$, returns the group $(S/I)^\times$ and a map $(S/I)^\times \to  S$ giving representatives. This is implemented only when $S$ is the maximal order.
intrinsic ResidueRingUnits(S::AlgEtQOrd,I::AlgEtQIdl) -> GrpAb,Map
{Given an order S and a fractional S-ideal, returns the group (S/I)^* and a map (S/I)^* -> S giving representatives. This is implemented only when S is the maximal order.}
    if not assigned I`ResidueRingUnits then
        require Order(I) eq S and I subset OneIdeal(S) : "I is not a proper S-ideal.";
        require IsMaximal(MultiplicatorRing(I)) : "the ideals needs to have maximal multiplicator ring.";
        maximal_order_case:=function(S,I)
        // if S is maximal
            vprintf ResidueRingUnits,2 : "ResidueRingUnits, maximal_order_case\n";
            test,I_asProd:=IsProductOfIdeals(I);
            assert test;
            A:=Algebra(S);
            n:=#I_asProd;
            _,embs:=Components(A);
            ray_res_rings:=[];
            ray_res_rings_maps:=[**];
            for i in [1..n] do
                vprintf ResidueRingUnits,2 : "ResidueRingUnits, maximal_order_case: %o-th/%o component...",i,n;
                IL:=I_asProd[i];
                OL:=Order(IL);
                assert IsMaximal(OL);
                vprintf ResidueRingUnits,2 : "...RayResidueRing...";
                R,r:=RayResidueRing(IL);
                vprintf ResidueRingUnits,2 : "...done\n";
                assert2 forall{ i : i,j in [1..Ngens(R)] | r(R.i+R.j) - r(R.i)*r(R.j) in IL };
                assert2 forall{ i : i in [1..Ngens(R)] | r(R.i) in OL };
                Append(~ray_res_rings,R);
                Append(~ray_res_rings_maps,r);
            end for;
            D,mRD,mDR:=DirectSum(ray_res_rings);
            vprintf ResidueRingUnits,2 : "ResidueRingUnits, maximal_order_case: DirectSum done\n";

            map_ResRing_S:=function(x)
                return &+[embs[i](ray_res_rings_maps[i](mDR[i](x))) : i in [1..n]];
            end function;

            map_S_ResRing:=function(y)
                comp:=Components(y);
                assert #ray_res_rings_maps eq #comp;
                return &+[mRD[i](comp[i]@@ray_res_rings_maps[i]) : i in [1..n]];
            end function;

            map:=map<D -> A | x:->map_ResRing_S(x) , y:->map_S_ResRing(y) >;
            vprintf ResidueRingUnits,2 : "ResidueRingUnits, maximal_order_case: map done\n";
            assert2 forall{ gen : gen in Generators(D) | (map(gen))@@map eq gen };
            assert2 forall{ gen : gen in Generators(D) | (map(gen)) in S };
            assert2 forall{ i : i,j in [1..Ngens(D)] | map(D.i+D.j) - map(D.i)*map(D.j) in I};
            return D,map;
        end function;

        if IsMaximal(S) then
            D,map:=maximal_order_case(S,I);
            I`ResidueRingUnits:=<D,map>;
        else
            OO:=MaximalOrder(Algebra(S));
            IOO:=OO!!I;
            O_I,map:=maximal_order_case(OO,IOO);
            gens:=ResidueRingUnitsSubgroupGenerators(I);
            assert2 forall{ g : g in gens | map(g@@map)-g in IOO };
            //assert2 forall{ g : g in gens | map(g@@map)-g in I };
            gens:=[ g@@map : g in gens ];
            D:=sub<O_I | gens >;
            I`ResidueRingUnits:=<D,map>;
        end if;
        assert2 forall{ gen : gen in Generators(D) | (map(gen))@@map eq gen };
        assert2 forall{ i : i,j in [1..Ngens(D)] | map(D.i+D.j) - map(D.i)*map(D.j) in I};
    end if;
    return Explode(I`ResidueRingUnits);
end intrinsic;

///ditto
intrinsic ResidueRingUnits(I::AlgEtQIdl) -> GrpAb,Map
{Returns the group (S/I)^* and a map (S/I)^* -> S, where S=Order(I) and the multiplicator ring of I is maximal.}
    return ResidueRingUnits(Order(I),I);
end intrinsic;

/// Given a fractional $S$-ideal $F$, returns generators of $(S/F)^\times$.
intrinsic ResidueRingUnitsSubgroupGenerators(F::AlgEtQIdl) -> SeqEnum[AlgEtQElt]
{Given a fractional S-ideal F, returns generators of (S/F)^*.}
    if not assigned F`residue_class_ring_unit_subgroup_generators then
        S:=Order(F);
        A:=Algebra(S);
        O:=MaximalOrder(A);
        Fm:=O!!F;
        l:=Factorization(Fm);
        l2:=[ <(S!!x[1]) meet S,x[2]>: x in l];
        primes:={x[1]:x in l2};
        primes:=[<x, Maximum([y[2]: y in l2 |y[1] eq x])>:x in primes];
        primes_powers:=[ p[1]^p[2] : p in primes ];
        elts:={};
        for i->a in primes do
            a1a2:=primes_powers[i];
            idp:=a[1];
            if #primes gt 1 then
                // coprime ideals, so we can use meet instead of *, which is faster
                rest:=&meet[ primes_powers[j] : j in [1..#primes] | j ne i];
            else
                rest:=OneIdeal(S);
            end if;
            //Compute primitive elt for residue field
            c:=PrimitiveElementResidueField(idp);
            e1:=ChineseRemainderTheorem(a1a2,rest,One(A),Zero(A));
            e2:=ChineseRemainderTheorem(a1a2,rest,Zero(A),One(A));
            c:=c*e1+e2;
            // the previous is equivalent to the following, but faster since we can 'reuse'e1 and e2 later.
            //c:=ChineseRemainderTheorem(a1a2,rest,c,One(A));
            assert2 c - ChineseRemainderTheorem(a1a2,rest,c,One(A)) in a1a2*rest;
            assert2 c - One(A) in rest;
            Include(~elts,c);
            b:=1;
            while b lt a[2] do
                M:=ZBasis(idp);
                M:=[1+x:x in M];
                for elt in M do
                    c:=elt*e1+e2;
                    // equivalent to:
                    //c:=ChineseRemainderTheorem((a1a2),rest,elt,One(A));
                    assert2 c - ChineseRemainderTheorem(a1a2,rest,elt,One(A)) in a1a2*rest;
                    Include(~elts,c);
                end for;
                b:=b*2;
                idp:=idp^2;
            end while;
        end for;
        assert2 forall{x : x in elts | x in S and not x in F};
        F`residue_class_ring_unit_subgroup_generators:=elts;
        vprintf AlgEtQPicardGroup, 2 :"residue_class_ring_unit_subgroup_generators:\n
                                         elts = %o\n",PrintSeqAlgEtQElt(Setseq(elts));
    end if;
    return F`residue_class_ring_unit_subgroup_generators;
end intrinsic;

IsPrincipal_prod_internal_old:=function( I , GRH )
//returns if the argument is a principal ideal; if so the function returns also the generator. It works only for products of ideals
    assert IsMaximal(Order(I)); //this function should be called only for ideals of the maximal order
    if #Generators(I) eq 1 then return true,Generators(I)[1]; end if;
    test,I_asProd:=IsProductOfIdeals(I);
    assert test; //this function should be called only for ideals of the maximal order, hence I is a product
    S:=Order(I);
    A:=Algebra(S);
    nf,embs:=Components(A);
    gen:=Zero(A);
    // the optional parameter Proof for the ClassGroup call, which by default is "Full"
    proof := GRH select "GRH" else "Full";
    for i in [1..#I_asProd] do
        IL:=I_asProd[i];
        OrdIL:=Order(IL);
        assert IsMaximal(OrdIL);
        //The next call is to prevent a bug of the in-built function IsPrincipal (which might have been corrected by now...).
        OL,oL:=ClassGroup(OrdIL : Proof := proof);
        testL,genL:=IsPrincipal(IL);
        assert2 (Zero(OL) eq (IL@@oL)) eq testL;
        if not testL then
            return false,_;
        end if;
        gen:=gen+embs[i](nf[i] ! genL);
    end for;
    assert2 gen*S eq I;
    I`Generators:=[gen];
    return true,gen;
end function;

IsPrincipal_prod_internal:=function( II , GRH )
//returns if the argument is a principal ideal; if so the function returns also the generator. It works only for products of ideals
    assert IsMaximal(Order(II)); //this function should be called only for ideals of the maximal order
    if #Generators(II) eq 1 then return true,Generators(II)[1]; end if;
    // we try with LLL
    elt:=ShortElement(II);
    if elt*Order(II) eq II then
            II`Generators:=[elt];
            return true,Generators(II)[1];
    end if;


    I,a:=SmallRepresentative(II); //a*II=I
    vprintf AlgEtQPicardGroup, 3:"IsPrincipal_prod_internal:\n
                                    ZBasis(II) = %o\n,ZBasis(I) = %o\n",
                                    PrintSeqAlgEtQElt(ZBasis(II)),PrintSeqAlgEtQElt(ZBasis(I));
    test,I_asProd:=IsProductOfIdeals(I);
    assert test; //this function should be called only for ideals of the maximal order, hence I is a product
    S:=Order(I);
    A:=Algebra(S);
    nf,embs:=Components(A);
    gen:=Zero(A);
    // the optional parameter Proof for the ClassGroup call, which by default is "Full"
    proof := GRH select "GRH" else "Full";
    for i in [1..#I_asProd] do
        IL:=I_asProd[i];
        OrdIL:=Order(IL);
        assert IsMaximal(OrdIL);
        //The next call is to prevent a bug of the in-built function IsPrincipal (which might have been corrected by now...).
        //Also if one wants to use the GRH bound rather than one needs to precompute the class groups, since IsPrincipal does not accepts varargs )
        OL,oL:=ClassGroup(OrdIL : Proof := proof);
        // the next line seems to be the most consuming part of the package.
        // It is 10/20 times slower than the RngOrdIdl version. Need to investigate
        testL,genL:=IsPrincipal(IL);
        assert2 (Zero(OL) eq (IL@@oL)) eq testL;
        if not testL then
            return false,_;
        end if;
        gen:=gen+embs[i](nf[i] ! genL);
    end for;
    assert2 gen*S eq I;

    gen:=gen/a;
    vprintf AlgEtQPicardGroup, 3:"IsPrincipal_prod_internal:\n
                            [gen,a] = %o\n",PrintSeqAlgEtQElt([gen,a]);
    II`Generators:=[gen];
    return true,gen;
end function;

/// Returns whether the given fractional ideals is a principal ideal (over its order of definition), and, if so, also a generator.  The parameter `GRH` (default false) can be set to true determines if the test is done under GRH.
intrinsic IsPrincipal(I1::AlgEtQIdl : GRH:=false )->BoolElt, AlgEtQElt
{Returns whethere the given fractional ideals is a principal ideal (over its order of definition), and, if so, also a generator.  The parameter GRH (default false) can be set to true determines if the test is done under GRH.}
    if not IsInvertible(I1) then return false,_; end if;
    if #Generators(I1) eq 1 then return true,Generators(I1)[1]; end if;
    // we try with LLL
    elt:=ShortElement(I1);
    if elt*Order(I1) eq I1 then
            I1`Generators:=[elt];
            return true,Generators(I1)[1];
    end if;

    S:=Order(I1);
    if IsMaximal(S) then
        return IsPrincipal_prod_internal(I1,GRH);
    end if;
    A:=Algebra(S);
    O:=MaximalOrder(A);
    F:=Conductor(S);
    FO:=O!!F;
    cop,I:=CoprimeRepresentative(I1,F);
    IO:=O!!I;
    is_princ_IO,gen_IO:=IsPrincipal_prod_internal(IO,GRH);
    if not is_princ_IO then
        return false,_;
    end if;
    R,r:=ResidueRingUnits(O,FO);
    if Order(R) eq 1 then
        assert2 gen_IO*S eq I;
        return true, gen_IO*cop^-1;
    end if;
    UO,uO:=UnitGroup(O : GRH:=GRH );
    Sgens:=ResidueRingUnitsSubgroupGenerators(F);
    B,b:=quo<R|[gen@@r : gen in Sgens]>;
    //B,b:=ResidueRingUnits(F); // it doesn't see that it is related to R.
    gens_UO_inB:=[ b(uO(UO.i)@@r) : i in [1..#Generators(UO)]  ];
    h:=hom<UO -> B | gens_UO_inB >;
    hUO:=Image(h);
    if not b(gen_IO@@r) in hUO then
        return false,_;
    end if;
    //now we know that IO is principal. let's find a generator
    UQ,qQ:=quo<UO|Kernel(h)>;  //UQ = O*/S*
    alpha:=hom<UQ -> B | [UQ.i@@qQ @uO @@r @b : i in [1..#Generators(UQ)]]>;
    is_princ,elt:=HasPreimage(gen_IO@@r@b,alpha);
    if is_princ then
        gen_I:=gen_IO*(A ! (elt@@qQ@uO))^-1;
        gen_I1:=gen_I*cop^-1;
        assert2 gen_I1*S eq I1;
        I1`Generators:=[gen_I1];
        return true,gen_I1;
    else
        return false, _;
    end if;
end intrinsic;

PicardGroup_prod_internal:=function( O , GRH )
//computes the PicardGroup of a product of order in a product of number fields and returns the group (as a direct product) and a sequence of representatives
    if assigned O`PicardGroup then
        return Explode(O`PicardGroup);
    end if;
    A:=Algebra(O);
    assert IsMaximal(O); // this function should be used only for maximal orders
    test,O_asProd:=IsProductOfOrders(O);
    assert test; //O must be a product of orders
    nf,embs:=Components(A);
    assert #nf eq #O_asProd;
    groups_maps_fields_maps:=[**];
    // the optional parameter Proof for the ClassGroup call, which by default is "Full"
    proof := GRH select "GRH" else "Full";
    for i in [1..#O_asProd] do
        OL:=O_asProd[i];
        GL,gL:=ClassGroup(OL : Proof := proof);
        assert2 forall{y : y in [gL(z) : z in GL] | MultiplicatorRing(y) eq OL};
        //this is a detector for bugs for the PicardGroup function.
        Append(~groups_maps_fields_maps,<GL,gL,nf[i],embs[i]>);
    end for;
    assert #groups_maps_fields_maps eq #nf;
    G,g,Gproj:=DirectSum([T[1] : T in groups_maps_fields_maps]);

    if #G eq 1 then
        from_G_to_ideals:=function(x)
            return OneIdeal(O);
        end function;
        from_ideals_to_G:=function(y)
            //the function is used only for max orders. no need to assert invertibility of y.
            assert Order(y) eq O;
            return Zero(G);
        end function;
        codomain:=Parent(OneIdeal(O));
        O`PicardGroup:=<G,map<G -> codomain | x:-> from_G_to_ideals(x) , y:->from_ideals_to_G(y) >>;
    else
        zerosinO:=[ Ideal(O,[ T[4](y) : y in Basis(T[2](Zero(T[1])),T[3])]) : T in groups_maps_fields_maps];
        assert2 &+zerosinO eq OneIdeal(O);
        geninO:=[]; //this will contain the the ideals of O corresponding to the generators of G
        for i in [1..#Generators(G)] do
            gen:=G.i;
            gens_inA:=[];
            for i in [1..#groups_maps_fields_maps] do
                T:=groups_maps_fields_maps[i];
                gLi:=T[2];
                idLi:=gLi(Gproj[i](gen));
                gens_inA:=gens_inA cat[T[4](g) : g in Basis(idLi,T[3])];
            end for;
            gen_O:=Ideal(O,gens_inA);
            // make the rep small
            gen_O:=SmallRepresentative(gen_O); //the ZBasis is already LLL-ed
            assert2 IsInvertible(gen_O); //test using colon ideal
            gen_O`IsInvertible:=true;
            gen_O`MultiplicatorRing:=O;
            // invertible, so 2-generated.
            TwoGeneratingSet(gen_O);
            Append(~geninO,gen_O);
        end for;
        assert #geninO eq #Generators(G);
        rep_idinA:= function(x)
            vprint AlgEtQPicardGroup, 2: "PicardGroup_prod_internal: rep_idinA\n";
            coeff:=Eltseq(x);
            id:=&*[geninO[i]^coeff[i] : i in [1..#coeff]];
            return id;
        end function;
        inverse_map:=function(id)
            vprint AlgEtQPicardGroup, 2: "PicardGroup_prod_internal: inverse_map\n";
            assert IsInvertible(id);
            if not IsIntegral(id) then
                id:=MakeIntegral(id);
            end if;
            test,id_asprod:=IsProductOfIdeals(id);
            assert test;
            return &+[g[i](id_asprod[i]@@groups_maps_fields_maps[i,2]) : i in [1..#id_asprod]];
        end function;
        Codomain:=Parent(rep_idinA(Zero(G)));
        mapGtoO:=map<G -> Codomain | rep:-> rep_idinA(rep) , y:->inverse_map(y) >;
        assert2 forall{a : a in Generators(G)| (mapGtoO(a))@@mapGtoO eq a};
        O`PicardGroup:=<G,mapGtoO>;
        vprintf AlgEtQPicardGroup, 2:"PicardGroup_prod_internal:\n
                              geninO = %o\n",[PrintSeqAlgEtQElt(ZBasis(I)) : I in geninO];
        O`PicardGroup:=<G,mapGtoO>;
    end if;
    return Explode(O`PicardGroup);
end function;

/// Returns the Picard group $\mathrm{Pic}(S)$ of the order $S$, which is not required to be maximal, and a map from the group to a set of representatives of the ideal classes.  The parameter `GRH` (default false) can be set to true determines if the test is done under GRH.
intrinsic PicardGroup( S::AlgEtQOrd : GRH:=false ) -> GrpAb, Map
{Returns the Picard group of the order S, which is not required to be maximal, and a map from the group to a set of representatives of the ideal classes.  The parameter GRH (default false) can be set to true determines if the test is done under GRH.}
    if assigned S`PicardGroup then
        return Explode(S`PicardGroup);
    end if;
    if IsMaximal(S) then
        return PicardGroup_prod_internal(S,GRH);
    end if;
    A:=Algebra(S);
    O:=MaximalOrder(A);
    GO,gO:=PicardGroup_prod_internal(O,GRH); //C, mC
    F:=Conductor(S);
    FO:=O!!F; //Fm
    gens_GO_in_S:=[]; //coprime with FO, in S and then meet S
    gens_GO_in_O:=[]; //coprime with FO, in O, Cgen
    if #GO gt 1 then
        for i in [1..#Generators(GO)] do
            I:=gO(GO.i); //already made small in the internal function
            c,cI:=CoprimeRepresentative(I,FO);
            //c,cI:=CoprimeRepresentative(I,MinimalInteger(FO)*O);
            cISmeetS:=(S!!cI) meet S;
            assert2 IsInvertible(cISmeetS); //test using colon ideal
            cISmeetS`IsInvertible:=true;
            cISmeetS`MultiplicatorRing:=S;
            // invertible, so 2-generated.
            TwoGeneratingSet(cISmeetS);
            Append(~gens_GO_in_S,cISmeetS);
            Append(~gens_GO_in_O,cI);//used in building relDglue
        end for;

        mGO_to_S:=function(rep)
            vprint AlgEtQPicardGroup, 3: "PicardGroup: mGO_to_S\n";
            coeff:=Eltseq(rep);
            idS:=&*[(gens_GO_in_S[i])^coeff[i] : i in [1..#coeff] ];
            return idS;
        end function;
        mGO_to_O:=function(rep)
            vprint AlgEtQPicardGroup, 3: "PicardGroup: mGO_to_O\n";
            coeff:=Eltseq(rep);
            assert #coeff eq #gens_GO_in_O;
            idO:=&*[(gens_GO_in_O[i])^coeff[i] : i in [1..#coeff] ];
            return idO;
        end function;
    else
        GO:=FreeAbelianGroup(0);
        mGO_to_S:=function(rep)
            idS:=OneIdeal(S);
            return idS;
        end function;
        mGO_to_O:=function(rep)
            idO:=OneIdeal(O);
            return idO;
        end function;
    end if;
    vprintf AlgEtQPicardGroup, 2:"PicardGroup:\n
                        gens_GO_in_S = %o\n
                        gens_GO_in_O = %o\n",
                        [PrintSeqAlgEtQElt(ZBasis(I)) : I in gens_GO_in_S],
                        [PrintSeqAlgEtQElt(ZBasis(I)) : I in gens_GO_in_O];


    R,r:=ResidueRingUnits(O,FO); // G, mG
    Sgens:=ResidueRingUnitsSubgroupGenerators(F); // ogens //generators in S of (S/F)*
    vprintf AlgEtQPicardGroup, 2: "PicardGroup: UO...";
    UO,uO:=UnitGroup(O : GRH:=GRH ); // Um, mUm
    vprintf AlgEtQPicardGroup, 2: "done\n";
    H:=FreeAbelianGroup(#Generators(GO));
    D, mRD, mHD, mDR, mDH := DirectSum(R,H); // D, mGD, mHD, mDG, mDH
    relDresidue:=[mRD(x@@r) : x in Sgens];
    relDunits:=[mRD(uO(x)@@r)  : x in Generators(UO)];
    // glue R and GO
    relDglue := [];
    assert #gens_GO_in_S eq #InvariantFactors(GO);
    for i in [1..#gens_GO_in_S] do
        is_princ, gen := IsPrincipal_prod_internal(gens_GO_in_O[i]^InvariantFactors(GO)[i],GRH);
        assert is_princ;
        Append(~relDglue,-mRD(gen@@r)+mHD(H.i*InvariantFactors(GO)[i]));
    end for;

    P, mDP := quo<D|relDresidue cat relDunits cat relDglue>;
    vprintf AlgEtQPicardGroup, 2: "PicardGroup: D = %o\n",D;
    vprintf AlgEtQPicardGroup, 2: "PicardGroup: P = %o\n",P;
    vprintf AlgEtQPicardGroup, 2: "PicardGroup: ker(mDP) = %o\n",Kernel(mDP);
    gens_P_in_D:=[P.i@@mDP : i in [1..#Generators(P)]];
    // Issue: the elements in gens_P_in_D might have very large coefficients, resulting in generators with
    // horrible Z-bases. The next few lines, minimize the coefficients that appear in front of the free generators.
        gens_P_in_D_small:=[];
        ker:=Kernel(mDP);
        ker_gens_inD:=[D!ker.i:i in [1..Ngens(ker)]];
        ker_free_gens_inD:=[g:g in ker_gens_inD|Order(g) eq 0];
        if #ker_free_gens_inD gt 0 then
            basis_mat:=Matrix([Eltseq(mDH(g)):g in ker_free_gens_inD]);
            basis_mat_LLL,T:=LLL(basis_mat); // basis_mat_LLL = T*basis_mat
            Latt:=LatticeWithBasis(basis_mat_LLL);
            vprintf AlgEtQPicardGroup, 2: "PicardGroup: basis_mat = \n%o\n",basis_mat;
            vprintf AlgEtQPicardGroup, 2: "PicardGroup: basis_mat_LLL = \n%o\n",basis_mat_LLL;
            for g in gens_P_in_D do
                elt:=Vector(Eltseq(-mDH(g)));
                vprintf AlgEtQPicardGroup, 2: "PicardGroup: ClosestVector to %o...",elt;
                c:=ClosestVector(Latt,elt);
                vprintf AlgEtQPicardGroup, 2: "done\n";
                c:=Eltseq(Vector(Coordinates(c))*T);
                g_small:=g+&+[c[i]*ker_free_gens_inD[i]:i in [1..#c]];
                // replace coeffs c of finite order N with N-c if better size
                c_small:=Eltseq(g_small);
                orders_Di:=[Order(D.i) : i in [1..Ngens(D)]];
                N_c_small:=[c_small[i]-orders_Di[i]:i in [1..#c_small]];
                c_better:=[Order(D.i) ne 0 and (Abs(N_c_small[i]) lt Abs(c_small[i])) select N_c_small[i] else c_small[i] : i in [1..#c_small]];
                g_small:=&+[c_better[i]*D.i:i in[1..Ngens(D)]];
                assert g-g_small in ker;
                vprintf AlgEtQPicardGroup, 2: "PicardGroup: g,g_small = %o\t%o\n",g,g_small;
                if &+[Abs(c):c in Eltseq(g)] lt &+[Abs(c):c in Eltseq(g_small)] then
                    Append(~gens_P_in_D_small,g);
                else
                    Append(~gens_P_in_D_small,g_small);
                end if;
            end for;
            gens_P_in_D:=gens_P_in_D_small;
        end if;
    // End of optimization

    if #P gt 1 then
        generators_ideals:=[];
        for igen->gen in gens_P_in_D do
            vprintf AlgEtQPicardGroup, 2: "PicardGroup: gens_P_in_D[%o]i = %o...\n",igen,gen;
            vprintf AlgEtQPicardGroup, 2: "\t\tcomputing id1 ...";
            id1:=(S!!( O*(r(mDR(gen))) )) meet S;
            vprintf AlgEtQPicardGroup, 2: "done: it has ZBasis of length %o\n",#Sprint(ZBasis(id1));
            vprintf AlgEtQPicardGroup, 2: "\t\tcomputing id2 ...";
            id2:=mGO_to_S(mDH(gen));
            vprintf AlgEtQPicardGroup, 2: "done: it has ZBasis of length %o\n",#Sprint(ZBasis(id2));
            vprintf AlgEtQPicardGroup, 2: "\t\tcomputing id1*id2 ...";
            gen_inS:=id1*id2;
            vprintf AlgEtQPicardGroup, 2: "done: it has ZBasis of length %o\n",#Sprint(ZBasis(gen_inS));
            // make the rep small
            vprintf AlgEtQPicardGroup, 2: "\t\tcalling SmallRepresentative ...";
            gen_inS:=SmallRepresentative(gen_inS); // the ZBasis is already LLL-ed
            vprintf AlgEtQPicardGroup, 2: "done\n";
            assert2 IsInvertible(gen_inS); //test using colon ideal
            gen_inS`IsInvertible:=true;
            gen_inS`MultiplicatorRing:=S;
            // invertible, so 2-generated.
            vprintf AlgEtQPicardGroup, 2: "\t\tcalling TwoGeneratingSet ...";
            TwoGeneratingSet(gen_inS);
            vprintf AlgEtQPicardGroup, 2: "done\n";
            Append(~generators_ideals,gen_inS);
        end for;
    else
        p:=map<P->Parent(OneIdeal(S)) | rep:->OneIdeal(S),
                                              id:->Zero(P)>;
        S`PicardGroup:=<P,p>;
        return Explode(S`PicardGroup);
    end if;

    vprint AlgEtQPicardGroup, 2: "PicardGroup: defining representative_picard_group ...";
    representative_picard_group := function(rep)
        vprint AlgEtQPicardGroup, 3: "PicardGroup: representative_picard_group\n";
        repseq := Eltseq(rep);
        return &*[generators_ideals[i]^repseq[i]:i in [1..#generators_ideals]];
    end function;
    vprint AlgEtQPicardGroup, 2: "done\n";


    vprint AlgEtQPicardGroup, 2: "PicardGroup: defining disc_log_picard_group ...";
    disc_log_picard_group:=function(id)
    // (crep*id)^-1 is coprime with F
        vprint AlgEtQPicardGroup, 3: "PicardGroup: disc_log_picard_group\n";
        crep:=1/(CoprimeRepresentative(id^-1,F));//here we check if id is invertible
        idO:=O!!(crep*id); //idO is coprime with FO
        GOrep:=idO@@gO;
        J:=mGO_to_O((H!Eltseq(GOrep))); //no minus signs, so J is coprime with FO
        assert2 IsCoprime(J,FO);
        prod:=(idO^-1)*J; //prod is corpime with FO...
        isprinc,elt:=IsPrincipal_prod_internal(prod,GRH);
        assert2 IsCoprime(elt*O,FO); //..hence elt is in the image r:R->Pic(S)
        Rrep:=elt@@r;
        rep_P:=mDP(-mRD(Rrep)+mHD(H!Eltseq(GOrep)));//[I]=-[xO meet S]+GOrep
        return rep_P;
    end function;
    vprint AlgEtQPicardGroup, 2: "done\n";

    vprint AlgEtQPicardGroup, 2: "PicardGroup: computing cod ...";
    cod:=Parent(representative_picard_group(Zero(P)));
    vprint AlgEtQPicardGroup, 2: "done\n";

    vprint AlgEtQPicardGroup, 2: "PicardGroup: defining the map p ...";
    p:=map<P -> cod | rep:->representative_picard_group(rep),
                       id:->disc_log_picard_group(id) >;
    vprint AlgEtQPicardGroup, 2: "done\n";
    S`PicardGroup:=<P,p>;
    vprintf AlgEtQPicardGroup, 2:"PicardGroup: computation is concluded\n
                        generators_ideals = %o\n",
                        [PrintSeqAlgEtQElt(ZBasis(I)) : I in generators_ideals];
    return P,p;
end intrinsic;

/// Given orders $S$ and $T$ such that $S \subseteq T$, returns the surjective extension map $\mathrm{Pic}(S) \to \mathrm{Pic}(T)$. The parameter `GRH`, default false, is passed to `PicardGroup`.
intrinsic ExtensionHomPicardGroups(S::AlgEtQOrd,T::AlgEtQOrd : GRH:=false)->Map
{Given orders S and T such that S subseteq T, returns the surjective extension map from PicardGroup(S) to PicardGroup(T). The vararg GRH, default false, is passed to PicardGroup.}
    require S subset T : "The second order is not an overorder of the first.";
    PS,pS:=PicardGroup(S:GRH:=GRH);
    PT,pT:=PicardGroup(T:GRH:=GRH);
    if S eq T then
        return IdentityHomomorphism(PS);
    else
        return hom<PS->PT|[(T!!(PS.i@pS))@@pT : i in [1..Ngens(PS)]]>;
    end if;
end intrinsic;

UnitGroup_prod_internal:=function(O, GRH)
    //returns the UnitGroup of a order which is a produc of orders
    if assigned O`UnitGroup then
        return Explode(O`UnitGroup);
    end if;
    assert IsMaximal(O); //this function should be used only for maximal orders
    test,O_asProd:=IsProductOfOrders(O);
    assert test; //the order must be a product
    A:=Algebra(O);
    nf,embs:=Components(A);
    idemA:=OrthogonalIdempotents(A);
    U_asProd:=[];
    u_asProd:=[**];
    for OL in O_asProd do
        U,u:=UnitGroup(OL : GRH:=GRH );
        Append(~U_asProd,U);
        Append(~u_asProd,u);
    end for;
    Udp,udp,proj_Udp:=DirectSum(U_asProd);
    gensinA:=[&+[embs[j](u_asProd[j](proj_Udp[j](Udp.i))) : j in [1..#U_asProd]] : i in [1..#Generators(Udp)] ];

    rep_inA:=function(rep)
        vprint AlgEtQPicardGroup, 2: "UnitGroup_prod_internal: rep_inA\n";
        coeff:=Eltseq(rep);
        return &*[gensinA[i]^coeff[i] : i in [1..#coeff]];
    end function;

    disc_log:=function(x)
        vprint AlgEtQPicardGroup, 2: "UnitGroup_prod_internal: disc_log\n";
        comp_x:=Components(A ! x);
        x_in_Udp:=&*[ udp[i](comp_x[i]@@u_asProd[i]) : i in [1..#comp_x] ];
        return x_in_Udp;
    end function;

    maptoA:=map<Udp -> A | rep :-> rep_inA(rep) , y :-> disc_log(y) >;
    O`UnitGroup:=<Udp,maptoA>;
    vprintf AlgEtQPicardGroup, 2:"UnitGroup_prod_internal:\n
                                    gensinA = %o\n",PrintSeqAlgEtQElt(gensinA);
    return Udp,maptoA;
end function;

/// Returns the unit group of the given order. The parameter `GRH` (default false) can be set to true determines if the test is done under GRH.
intrinsic UnitGroup(S::AlgEtQOrd : GRH:=false ) -> GrpAb, Map
{Returns the unit group of the given order. The parameter GRH (default false) can be set to true determines if the test is done under GRH.}
    if assigned S`UnitGroup then
        return Explode(S`UnitGroup);
    end if;
    if IsMaximal(S) then return UnitGroup_prod_internal(S,GRH); end if;
    A:=Algebra(S);
    require assigned A`Components: "the order must lie in a product of number fields";
    O:=MaximalOrder(Algebra(S));
    UO,uO:=UnitGroup_prod_internal(O, GRH);
    F:=Conductor(S);
    FO:=O!!F;

    //Let's build B=(O/FO)^*/(S/F)^*
    R,r:=ResidueRingUnits(O,FO);
    gens_SF:=ResidueRingUnitsSubgroupGenerators(F);
    B,b:=quo<R| [ a@@r : a in gens_SF ]>;

    img_gensUO_in_B:=[ b(uO(UO.i)@@r) : i in [1..#Generators(UO)] ];
    m:=hom<UO -> B | img_gensUO_in_B >;
    P:=Kernel(m);
    gens_P_in_A:=[uO(UO ! P.i) : i in [1..#Generators(P)] ];
    //p_codomain:=Parent(gens_P_in_A[1]);
    p_codomain:=A;

    map_P_to_S:=function(rep)
        vprint AlgEtQPicardGroup, 2: "UnitGroup: map_P_to_S\n";
        coeff:=Eltseq(rep);
        assert #coeff eq #gens_P_in_A;
        elt:=&*[gens_P_in_A[i]^coeff[i] : i in [1..#coeff]];
        assert2 elt in S and elt^-1 in S;
        return elt;
    end function;

    map_S_to_P:=function(y)
        vprint AlgEtQPicardGroup, 2: "UnitGroup: map_S_to_P\n";
        assert2 y in S and y^-1 in S;
        assert2 m(y@@uO) eq Zero(B);
        elt := P ! (y@@uO);
        return elt;
    end function;
    p:=map<P -> p_codomain | x:->map_P_to_S(x), y:->map_S_to_P(y)  >;
    S`UnitGroup:=<P,p>;
    vprintf AlgEtQPicardGroup, 2:"UnitGroup:\n
                                gens_P_in_A = %o\n",PrintSeqAlgEtQElt(gens_P_in_A);
    return P,p;
end intrinsic;

/// Given two fractional ideals $I$ and $J$ over the same order, checks whether $I=x\cdot J$, for some unit $x$. If so, $x$ is also returned. The parameter `GRH` (default false) can be set to true determines if the test is done under GRH.
intrinsic IsIsomorphic(I::AlgEtQIdl, J::AlgEtQIdl : GRH:=false ) -> BoolElt, AlgEtQElt
{Given two fractional ideals I and J over the same order, checks whether I=x*J, for some unit x. If so, x is also returned. The parameter GRH (default false) can be set to true determines if the test is done under GRH.}
    test:=IsWeakEquivalent(I,J); //if so I=(I:J)*J and (I:J) is invertible in its MultiplicatorRing
    if test then
        S:=MultiplicatorRing(I);
        IS:=S!!I;
        JS:=S!!J;
        CIJ:=ColonIdeal(IS,JS);
        test2,x:= IsPrincipal(CIJ : GRH:=GRH );
        if test2 then
            return test2,x;
        assert2 I eq x*J;
        else
            return false, _ ;
        end if;
    else
        return false , _ ;
    end if;
end intrinsic;

/* TESTS

    printf "### Testing PicardGroup and UnitGroup:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    SetAssertions(2);
    _<x>:=PolynomialRing(Integers());

    // very fast
    f:=(x^4+16);
    A:=EtaleAlgebra(f);
    E:=EquationOrder(A);
    oo:=FindOverOrders(E);
    for S in oo do
        P,p:=PicardGroup(S);
        U,u:=UnitGroup(S);
        printf ".";
    end for;

    //AttachSpec("~/packages_github/AlgEt/spec");
    _<x>:=PolynomialRing(Integers());
    f:=x^4-1000*x^3-1000*x^2-1000*x-1000;
    SetAssertions(1);
    SetVerbose("AlgEtQPicardGroup",1);
    SetVerbose("AlgEtQIdl",1);
    SetVerbose("ShortEltSmallRep",1);
    A:=EtaleAlgebra(f);
    E:=EquationOrder(A);
    P,p:=PicardGroup(E : GRH:=true);
    U,u:=UnitGroup(E : GRH:=true);
    assert #P eq 3548000;
    printf ".";
    SetAssertions(1);

    //AttachSpec("~/packages_github/AlgEt/spec");
    _<x>:=PolynomialRing(Integers());
    f:=x^8+16;
    A:=EtaleAlgebra(f);
    OA:=MaximalOrder(A);
    P:=PrimesAbove(2*OA)[1];
    E:=Order(ZBasis(P^2));
    assert not IsMaximal(E);
    P:=Conductor(E);
    assert IsPrime(P);
    assert #ResidueRingUnits(P) eq Index(E,P)-1;
    assert forall{ i : i in [1..20] | #ResidueRingUnits(Pi) eq Index(E,Pi) - Index(P,Pi) where Pi:=P^i};

    //AttachSpec("~/AlgEt/spec");
    _<x>:=PolynomialRing(Integers());
    f:=x^4-1000*x^3-1000*x^2-1000*x-1000;
    SetAssertions(1);
    SetVerbose("AlgEtQPicardGroup",1);
    SetVerbose("AlgEtQIdl",1);
    SetVerbose("ShortEltSmallRep",1);
    A:=EtaleAlgebra(f);
    E:=EquationOrder(A);
    P,p:=PicardGroup(E : GRH:=true);
    O:=MaximalOrder(A);
    C,c:=PicardGroup(O);
    ext:=ExtensionHomPicardGroups(E,O);
    for i in [1..100] do
        printf ".";
        x:=Random(P); y:=Random(P);
        assert ext(x+y) eq ext(x)+ext(y);
    end for;
    printf " all good!";
*/
