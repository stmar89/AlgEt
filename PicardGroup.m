/* vim: set syntax=magma :*/

//freeze;

/////////////////////////////////////////////////////
// Picard Group of orders in etale algebras over \Q
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose AlgEtPicardGroup, 3;

declare attributes AlgEtOrd:PicardGroup,
                            UnitGroup;
declare attributes AlgEtIdl:residue_class_ring_unit_subgroup_generator;

intrinsic ResidueRingUnits(S::AlgEtOrd,I::AlgEtIdl) -> GrpAb,Map
{returns the group (S/I)^* and a map (S/I)^* -> S. It is required S to be maximal }
    //the following code works only for maximal orders in etale algebras
    require IsMaximal(S): "implemented only for the maximal order";
    test,I_asProd:=IsProductOfIdeals(I);
    assert test;
    A:=Algebra(S);
    n:=#I_asProd;
    _,embs:=Components(A);
    ray_res_rings:=[];
    ray_res_rings_maps:=[**];
    for i in [1..n] do
        IL:=I_asProd[i];
        OL:=Order(IL);
        assert IsMaximal(OL);
        R,r:=RayResidueRing(IL);
        Append(~ray_res_rings,R);
        Append(~ray_res_rings_maps,r);
    end for;
    D,mRD,mDR:=DirectSum(ray_res_rings);

    map_ResRing_S:=function(x)
        return &+[embs[i](ray_res_rings_maps[i](mDR[i](x))) : i in [1..n]];
    end function;

    map_S_ResRing:=function(y)
        comp:=Components(y);
        assert #ray_res_rings_maps eq #comp;
        return &+[mRD[i](comp[i]@@ray_res_rings_maps[i]) : i in [1..n]];
    end function;

    map:=map<D -> A | x:->map_ResRing_S(x) , y:->map_S_ResRing(y) >;
    assert2 forall{ gen: gen in Generators(D) | (map(gen))@@map eq gen };
    return D,map;
end intrinsic;

// USE INSTEAD PrimitiveElementResidueField from Quotients.m
// residue_class_field_primitive_element := function(P)
// //given a maximal ideal P in S, returns a generator of (S/P)* as an element of S;
//     S:=Order(P);
//     Q,m:=ResidueRing(S,P);
//     ord:=#Q-1; //ord = #(S/P)-1;
//     assert IsPrimePower(#Q); // #(S/P) must be a prime power
//     proper_divisors_ord:=Exclude(Divisors(ord),ord);
//     repeat
//         repeat 
//             a:=Random(Q);
//         until a ne Zero(Q);
//     until forall{f : f in proper_divisors_ord | m((a@@m)^f) ne m(One(Algebra(S)))};
//     assert2 (m((a@@m)^ord) eq m(One(Algebra(S)))); 
//     return a@@m;
// end function;

residue_class_ring_unit_subgroup_generators:=function(F)
// determine generators of the subgroup of (S/F)^* as elements of A=Algebra(S)
    if not assigned F`residue_class_ring_unit_subgroup_generator then
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
                rest:=&*[ primes_powers[j] : j in [1..#primes] | j ne i];
            else
                rest:=OneIdeal(S);
            end if;
            //Compute primitive elt for residue field
            //c:=residue_class_field_primitive_element(idp);
            c:=PrimitiveElementResidueField(idp);
            e1:=ChineseRemainderTheorem(a1a2,rest,One(A),Zero(A));
            e2:=ChineseRemainderTheorem(a1a2,rest,Zero(A),One(A));
            c:=c*e1+e2;
            //c:=ChineseRemainderTheorem(a1a2,rest,c,One(A));
            assert2 c - ChineseRemainderTheorem(a1a2,rest,c,One(A)) in a1a2*rest;
            Include(~elts,c);
            b:=1;
            while b lt a[2] do
                M:=ZBasis(idp);
                M:=[1+x:x in M];
                for elt in M do
                    c:=elt*e1+e2;
                    //c:=ChineseRemainderTheorem((a1a2),rest,elt,One(A));
                    assert2 c - ChineseRemainderTheorem(a1a2,rest,elt,One(A)) in a1a2*rest;
                    Include(~elts,c);
                end for;
                b:=b*2;
                idp:=idp^2;
            end while;
        end for;
        assert2 forall{x : x in elts | x in S and not x in F};
        F`residue_class_ring_unit_subgroup_generator:=elts;
        vprintf AlgEtPicardGroup, 1:"residue_class_ring_unit_subgroup_generator:\n
                                         elts = %o\n",PrintSeqAlgEtElt(Setseq(elts));
    end if;
	return F`residue_class_ring_unit_subgroup_generator ;
end function;

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
    if GRH then
        SetClassGroupBounds("GRH");
    end if;
    for i in [1..#I_asProd] do
        IL:=I_asProd[i];
        OrdIL:=Order(IL);
        assert IsMaximal(OrdIL);
        //The next call is to prevent a bug of the in-built function IsPrincipal (which might have been corrected by now...).
        //Also if one wants to use the GRH bound rather than one needs to precompute the class groups, since IsPrincipal does not accepts varargs )       
        OL,oL:=ClassGroup(OrdIL);
        // the next line seems to be the most time consuming part of the package. 
        // It is 10/20 times slower than the RngOrdIdl version. Need to investigate
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

    I,a:=SmallRepresentative(II); //a*II=I
    vprintf AlgEtPicardGroup, 2:"IsPrincipal_prod_internal:\n
                                    ZBasis(II) = %o\n,ZBasis(I) = %o\n",
                                    PrintSeqAlgEtElt(ZBasis(II)),PrintSeqAlgEtElt(ZBasis(I));
    test,I_asProd:=IsProductOfIdeals(I);
    assert test; //this function should be called only for ideals of the maximal order, hence I is a product
    S:=Order(I);
    A:=Algebra(S);
    nf,embs:=Components(A);
    gen:=Zero(A);
    if GRH then
        SetClassGroupBounds("GRH");
    end if;
    for i in [1..#I_asProd] do
        IL:=I_asProd[i];
        OrdIL:=Order(IL);
        assert IsMaximal(OrdIL);
        //The next call is to prevent a bug of the in-built function IsPrincipal (which might have been corrected by now...).
        //Also if one wants to use the GRH bound rather than one needs to precompute the class groups, since IsPrincipal does not accepts varargs )       
        OL,oL:=ClassGroup(OrdIL);
        // the next line seems to be the most time consuming part of the package. 
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
    vprintf AlgEtPicardGroup, 2:"IsPrincipal_prod_internal:\n
                            [gen,a] = %o\n",PrintSeqAlgEtElt([gen,a]); 
    II`Generators:=[gen];
    return true,gen;
end function;

intrinsic IsPrincipal(I1::AlgEtIdl : GRH:=false )->BoolElt, AlgAssElt
{   
    Return if the argument is a principal ideal; if so the function returns also the generator.
    The optional argument "GRH" decides wheter the bound for the IsPrincipal test should be conditional. The default value is "false".
}
//wouldn't an LLL test be faster?
    if not IsInvertible(I1) then return false,_; end if;
    if #Generators(I1) eq 1 then return true,Generators(I1)[1]; end if;

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
    Sgens:=residue_class_ring_unit_subgroup_generators(F);
    B,b:=quo<R|[gen@@r : gen in Sgens]>;
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
    if GRH then
        SetClassGroupBounds("GRH");
    end if;
    for i in [1..#O_asProd] do
        OL:=O_asProd[i];
        GL,gL:=ClassGroup(OL);
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
        return G,map<G -> codomain | x:-> from_G_to_ideals(x) , y:->from_ideals_to_G(y) >;
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
            gen_O:=SmallRepresentative(gen_O);
            assert2 IsInvertible(gen_O); //test using colon ideal
            gen_O`IsInvertible:=true;
            gen_O`MultiplicatorRing:=O;
            // invertible, so 2-generated.
            TwoGeneratingSet(gen_O);
            Append(~geninO,gen_O);
        end for;
        assert #geninO eq #Generators(G);      
        rep_idinA:= function(x)
            vprint AlgEtPicardGroup, 2: "PicardGroup_prod_internal: rep_idinA\n";
            coeff:=Eltseq(x);
            id:=&*[geninO[i]^coeff[i] : i in [1..#coeff]];
            return id;
        end function;
        inverse_map:=function(id)
            vprint AlgEtPicardGroup, 2: "PicardGroup_prod_internal: inverse_map\n";
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
        vprintf AlgEtPicardGroup, 2:"PicardGroup_prod_internal:\n
                              geninO = %o\n",[PrintSeqAlgEtElt(ZBasis(I)) : I in geninO];
        return G,mapGtoO;
    end if;
end function;

intrinsic PicardGroup( S::AlgEtOrd : GRH:=false ) -> GrpAb, Map
{
    Return the PicardGroup of the order S, which is not required to be maximal, and a map from the PicardGroup to a set of representatives of the ideal classes
    The optional argument "GRH" decides the bound for the computations of the ClassGroup and UnitGroup of the maximal order. The default value is "false".
}
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
            vprint AlgEtPicardGroup, 2: "PicardGroup: mGO_to_S\n";
            coeff:=Eltseq(rep);
            idS:=&*[(gens_GO_in_S[i])^coeff[i] : i in [1..#coeff] ];
            return idS;
        end function;
        mGO_to_O:=function(rep) 
            vprint AlgEtPicardGroup, 2: "PicardGroup: mGO_to_O\n"; 
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
    vprintf AlgEtPicardGroup, 2:"PicardGroup:\n
                        gens_GO_in_S = %o\n
                        gens_GO_in_O = %o\n",
                        [PrintSeqAlgEtElt(ZBasis(I)) : I in gens_GO_in_S],
                        [PrintSeqAlgEtElt(ZBasis(I)) : I in gens_GO_in_O];


    R,r:=ResidueRingUnits(O,FO); // G, mG
    Sgens:=residue_class_ring_unit_subgroup_generators(F); // ogens //generators in S of (S/F)*
    UO,uO:=UnitGroup(O : GRH:=GRH ); // Um, mUm 
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
    gens_P_in_D:=[P.i@@mDP : i in [1..#Generators(P)]];
    if #P gt 1 then
        generators_ideals:=[];
        for gen in gens_P_in_D do
            id1:=(S!!( O*(r(mDR(gen))) )) meet S;
            id2:=mGO_to_S(mDH(gen));
            gen_inS:=id1*id2;
            // make the rep small
            gen_inS:=SmallRepresentative(gen_inS);
            assert2 IsInvertible(gen_inS); //test using colon ideal
            gen_inS`IsInvertible:=true;
            gen_inS`MultiplicatorRing:=S;
            // invertible, so 2-generated.
            TwoGeneratingSet(gen_inS);
            Append(~generators_ideals,gen_inS);
        end for;
    else
        return P,map<P->Parent(OneIdeal(S)) | rep:->OneIdeal(S),
                                              id:->Zero(P)>;
    end if;

    representative_picard_group := function(rep)
        vprint AlgEtPicardGroup, 2: "PicardGroup: representative_picard_group\n"; 
        repseq := Eltseq(rep);
        return &*[generators_ideals[i]^repseq[i]:i in [1..#generators_ideals]];
    end function;


    disc_log_picard_group:=function(id)
    // (crep*id)^-1 is coprime with F
        vprint AlgEtPicardGroup, 2: "PicardGroup: disc_log_picard_group\n"; 
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

    cod:=Parent(representative_picard_group(Zero(P)));
    p:=map<P -> cod | rep:->representative_picard_group(rep),
                       id:->disc_log_picard_group(id) >;
    S`PicardGroup:=<P,p>;
    vprintf AlgEtPicardGroup, 2:"PicardGroup:\n
                        generators_ideals = %o\n",
                        [PrintSeqAlgEtElt(ZBasis(I)) : I in generators_ideals];
    return P,p;
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
        vprint AlgEtPicardGroup, 2: "UnitGroup_prod_internal: rep_inA\n"; 
        coeff:=Eltseq(rep);
        return &*[gensinA[i]^coeff[i] : i in [1..#coeff]];
    end function;

    disc_log:=function(x)
        vprint AlgEtPicardGroup, 2: "UnitGroup_prod_internal: disc_log\n"; 
        comp_x:=Components(A ! x);
        x_in_Udp:=&*[ udp[i](comp_x[i]@@u_asProd[i]) : i in [1..#comp_x] ];
        return x_in_Udp;
    end function;

    maptoA:=map<Udp -> A | rep :-> rep_inA(rep) , y :-> disc_log(y) >;
    O`UnitGroup:=<Udp,maptoA>;
    vprintf AlgEtPicardGroup, 2:"UnitGroup_prod_internal:\n
                                    gensinA = %o\n",PrintSeqAlgEtElt(gensinA);
    return Udp,maptoA;
end function;

intrinsic UnitGroup(S::AlgEtOrd : GRH:=false ) -> GrpAb, Map
{   
    Return the unit group of a order in a etale algebra
    The optional argument "GRH" decides the bound for the computation of the unit group of the maximal order. The default value is "false".
}
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
    gens_SF:=residue_class_ring_unit_subgroup_generators(F);
    B,b:=quo<R| [ a@@r : a in gens_SF ]>;

    img_gensUO_in_B:=[ b(uO(UO.i)@@r) : i in [1..#Generators(UO)] ];
    m:=hom<UO -> B | img_gensUO_in_B >;
    P:=Kernel(m);
    gens_P_in_A:=[uO(UO ! P.i) : i in [1..#Generators(P)] ];
    //p_codomain:=Parent(gens_P_in_A[1]);
    p_codomain:=A;

    map_P_to_S:=function(rep)
        vprint AlgEtPicardGroup, 2: "UnitGroup: map_P_to_S\n"; 
        coeff:=Eltseq(rep);
        assert #coeff eq #gens_P_in_A;
        elt:=&*[gens_P_in_A[i]^coeff[i] : i in [1..#coeff]];
        assert2 elt in S and elt^-1 in S;
        return elt;
    end function;

    map_S_to_P:=function(y)
        vprint AlgEtPicardGroup, 2: "UnitGroup: map_S_to_P\n"; 
        assert2 y in S and y^-1 in S;
        assert2 m(y@@uO) eq Zero(B);
        elt := P ! (y@@uO);
        return elt;
    end function;
    p:=map<P -> p_codomain | x:->map_P_to_S(x), y:->map_S_to_P(y)  >;
    S`UnitGroup:=<P,p>;
    vprintf AlgEtPicardGroup, 2:"UnitGroup:\n
                                gens_P_in_A = %o\n",PrintSeqAlgEtElt(gens_P_in_A);
    return P,p;
end intrinsic;

intrinsic IsIsomorphic(I::AlgEtIdl, J::AlgEtIdl : GRH:=false ) -> BoolElt, AlgAssElt
{
    Checks if I=x*J, for some x. If so, also x is returned
    The optional argument "GRH" decides wheter the bound for the IsPrincipal test should be conditional. The default value is "false"
}
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

/*TEST
	
	AttachSpec("~/packages_github/AlgEt/spec");
	SetAssertions(1);
	_<x>:=PolynomialRing(Integers());

	// very fast
    f:=(x^4+16);
	A:=EtaleAlgebra(f);
	E:=EquationOrder(A);
	oo:=FindOverOrders(E);
    t0:=Cputime();
    for S in oo do
        P,p:=PicardGroup(S);
        U,u:=UnitGroup(S);	
    end for;
    Cputime(t0);

    // with Profiler ~16 sec with GRH:=true
	AttachSpec("~/packages_github/AlgEt/spec");
	SetAssertions(1);
	_<x>:=PolynomialRing(Integers());
    	SetProfile(true);
    for i in [1..10] do
        f:=x^4-1000*x^3-1000*x^2-1000*x-1000;
        A:=EtaleAlgebra(f);
        E:=EquationOrder(A);
        t0:=Cputime();
            P,p:=PicardGroup(E : GRH:=true);
            //U,u:=UnitGroup(E : GRH:=true);	
        Cputime(t0);
    end for;
	SetProfile(false);
    ProfilePrintByTotalTime(ProfileGraph());

    quit; 
    git pull;
    sleep 3; 
    magma;
	AttachSpec("~/packages_github/AlgEt/spec");
    //SetOutputFile("")
	_<x>:=PolynomialRing(Integers());
    f:=x^4-1000*x^3-1000*x^2-1000*x-1000;
    SetClassGroupBounds("GRH");
    SetVerbose("AlgEtPicardGroup",1);
    SetVerbose("AlgEtIdl",2);
    SetVerbose("ShortEltSmallRep",1);
    for i in [1..10^1] do
        //"NF";
        //time P,p:=PicardGroup(EquationOrder(NumberField(f)));
        //assert #P eq 3548000;
        "Et";
        time P,p:=PicardGroup(EquationOrder(EtaleAlgebra(f)) : GRH:=true);
        assert #P eq 3548000;
        "\n";
    end for;

    AttachSpec("~/packages_github/AlgEt/spec");
    _<x>:=PolynomialRing(Integers());
    f:=x^6 - 3*x^5 - 3*x^4 + 65*x^3 - 48*x^2 - 768*x + 4096;
    A:=EtaleAlgebra(f);
    gensT:=[
        <[ 1, 0 ], [ 1/9, 5/6, 1, 41/18 ]>,
        <[ 0, 1 ], [ 0, 1, 0, 0 ]>,
        <[ 0, 0 ], [ 8/9, 11/6, 4/3, 133/18 ]>,
        <[ 0, 0 ], [ 0, 8/3, 7/3, 29/3 ]>,
        <[ 0, 0 ], [ 0, 0, 3, 3 ]>,
        <[ 0, 0 ], [ 0, 0, 0, 18 ]>
    ];
    gensT:=[ A ! g : g in gensT ];
    T:=Order(gensT);
    #PicardGroup(T); // this used to trigger a bug in CRT. now fixed



*/
