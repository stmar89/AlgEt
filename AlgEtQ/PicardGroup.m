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

// Who writes what:
//  - PicardGroup and UnitGroup: the intrinsics with the same names, as <group, map, GRH>, where GRH
//    says whether GRH was assumed (see valid_for). The helpers that compute them return values and
//    store nothing.
//  - conductor_sequence: the function conductor_sequence (see conductor_sequence_fmt).
import "Ord.m" : all_but_one;

declare attributes AlgEtQOrd: PicardGroup,
                              UnitGroup,
                              conductor_sequence;

// Who writes what:
//  - residue_class_ring_unit_subgroup_generators: ResidueRingUnitsSubgroupGenerators;
//  - ResidueRingUnits: the intrinsic ResidueRingUnits, as <group, map>;
//  - SmallPowers: small_power; small_prod_cache: small_prod_of_powers; compact_powers:
//    compact_power.
declare attributes AlgEtQIdl: residue_class_ring_unit_subgroup_generators,
                              ResidueRingUnits,
                              SmallPowers,
                              small_prod_cache,
                              compact_powers;

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
/// 1 \to R^\times \to \mathcal{O}_A^\times \to \dfrac{\left( \mathcal{O}_A/\mathfrak{f} \right)^\times}{\left( R/\mathfrak{f} \right)^\times} \to \mathrm{Pic}(R)\to \mathrm{Pic}(\mathcal{O}_A) \to 1,
/// ```
/// where $\mathfrak{f}$ is the conductor of $R$, the first, second and third map are the natural maps, the fourth is induced by $a \mapsto (a\mathcal{O}_A \cap R)$ and the last one is induced by the extension map $I \mapsto I\mathcal{O}_A$.
///
/// The following intrinsics allow to compute all the elements of the sequence.
/// If one wants to perform the computation under the GRH bound one can either use the appropriate parameters or use `SetClassGroupBounds("GRH")` as for number fields.



/// Given an order $S$ and a fractional $S$-ideal $I$, returns the group $(S/I)^\times$ and a map $(S/I)^\times \to  S$ giving representatives. This is implemented only when $S$ is the maximal order.
intrinsic ResidueRingUnits(S::AlgEtQOrd,I::AlgEtQIdl) -> GrpAb,Map
{Given an order S and a fractional S-ideal, returns the group (S/I)^* and a map (S/I)^* -> S giving representatives. This is implemented only when S is the maximal order.}
    if not assigned I`ResidueRingUnits then
        require Order(I) eq S and I subset OneIdeal(S) : "I is not a proper S-ideal.";
        require IsMaximal(MultiplicatorRing(I)) : "the ideal needs to have maximal multiplicator ring.";
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

/// Given an integral $S$-ideal $F$ which is also an ideal of the maximal order of the algebra (for example, the conductor of $S$), returns generators of $(S/F)^\times$.
intrinsic ResidueRingUnitsSubgroupGenerators(F::AlgEtQIdl) -> SeqEnum[AlgEtQElt]
{Given an integral S-ideal F which is also an ideal of the maximal order of the algebra (for example, the conductor of S), returns generators of (S/F)^*.}
    if not assigned F`residue_class_ring_unit_subgroup_generators then
        S:=Order(F);
        require F subset OneIdeal(S) and IsMaximal(MultiplicatorRing(F)) : "F must be an integral ideal of S which is also an ideal of the maximal order.";
        A:=Algebra(S);
        O:=MaximalOrder(A);
        Fm:=O!!F;
        l:=Factorization(Fm);
        // the primes of S below the primes of O in the factorization: PrimeBelow returns the stored
        // primes of S, already marked as prime, so their residue fields are computed once (and
        // IsPrime is not needed)
        l2:=[ <PrimeBelow(S,x[1]),x[2]>: x in l];
        primes:={x[1]:x in l2};
        primes:=[<x, Maximum([y[2]: y in l2 |y[1] eq x])>:x in primes];
        primes_powers:=[ p[1]^p[2] : p in primes ];
        // rests[i] is the intersection of the prime powers other than the i-th (they are coprime,
        // so meet is the product, and it is faster), from the prefixes and the suffixes: see
        // all_but_one
        rests:=all_but_one(primes_powers,func< I,J | I meet J >,OneIdeal(S));
        elts:={};
        for i->a in primes do
            a1a2:=primes_powers[i];
            idp:=a[1];
            rest:=rests[i];
            //Compute primitive elt for residue field
            c:=PrimitiveElementResidueField(idp);
            e1:=ChineseRemainderTheorem(a1a2,rest,One(A),Zero(A));
            // e2 = 1-e1 is 0 modulo a1a2 and 1 modulo rest: no second CRT is needed
            e2:=One(A)-e1;
            c:=c*e1+e2;
            // the previous is equivalent to the following, but faster since we can 'reuse' e1 and
            // e2 later.
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

class_group_internal:=function(OL, proof)
// ClassGroup of the maximal order OL of a number field, with the parameter Proof (proof = "GRH" or
// "Full"). The orders of quadratic fields (type RngQuad, e.g. MaximalOrder(QuadraticField(d))) do
// not accept the parameter: Magma computes their class groups with binary quadratic forms.
    if Type(OL) eq RngQuad then
        return ClassGroup(OL);
    end if;
    return ClassGroup(OL : Proof := proof);
end function;


//------------
// Helpers: caches and GRH, closest vectors, invertible ideals
//------------

valid_for:=function(cached_GRH, GRH)
// Whether a stored result, computed with the parameter GRH equal to cached_GRH, can be returned
// when GRH is asked: a result computed without assuming GRH holds also under GRH, while one
// computed under GRH is computed again when it is asked without GRH.
    return GRH or not cached_GRH;
end function;

closest_vector_data:=function(M)
// Given an integer matrix M with linearly independent rows, returns <L,T>, where L is the lattice
// spanned by the rows of M with the LLL-reduced basis T*M: see closest_vector_coefficients.
    Lmat,T:=LLL(M); // Lmat = T*M
    return <LatticeWithBasis(Lmat),T>;
end function;

closest_vector_coefficients:=function(lat, v)
// Given <L,T> as returned by closest_vector_data(M), and a vector v, returns the integer
// coefficients c such that c*M is a closest vector to v in the lattice spanned by the rows of M.
    return Eltseq(Vector(Coordinates(ClosestVector(lat[1],v)))*lat[2]);
end function;

two_generated_invertible:=function(I, S)
// Given an invertible fractional S-ideal I, marks it as invertible with multiplicator ring S,
// stores two generators of it (see TwoGeneratingSet), and returns it.
    assert2 IsInvertible(I); // test using colon ideal
    I`IsInvertible:=true;
    I`MultiplicatorRing:=S;
    TwoGeneratingSet(I);
    return I;
end function;

//------------
// Small representatives of power products
//------------

small_power:=function(I, e)
// Given an invertible ideal I and a non-zero integer e, returns an ideal in the same ideal class as
// I^e, with a small Z-basis: every intermediate ideal is replaced by a small representative, so
// that the Z-bases do not blow up. The output is stored in the associative array I`SmallPowers,
// together with the powers met along the way (I^(e-1) for e odd, I^(e div 2) for e even, ..., and a
// small representative of I^-1 with key -1 for e negative): the representatives of many ideal
// classes, as needed for example by ICM, use the same powers of the generators of the Picard group
// over and over. At most 2*n powers are stored if the exponents are smaller than n.
    if e eq 1 then
        return I;
    end if;
    if not assigned I`SmallPowers then
        I`SmallPowers:=AssociativeArray();
    end if;
    is_def,J:=IsDefined(I`SmallPowers,e);
    if not is_def then
        if e eq -1 then
            J:=SmallRepresentative(I^-1);
        elif e lt 0 then
            return $$($$(I,-1),-e); // stored in the small representative of I^-1
        elif IsOdd(e) then
            J:=SmallRepresentative($$(I,e-1)*I);
        else
            H:=$$(I,e div 2);
            J:=SmallRepresentative(H*H);
        end if;
        J`IsInvertible:=true;
        I`SmallPowers[e]:=J;
    end if;
    return J;
end function;

small_prod_of_powers:=function(S, gens, exps)
// Given invertible S-ideals gens and integers exps, returns an ideal in the same ideal class as
// &*[gens[i]^exps[i]], with a small Z-basis: the powers are computed by small_power, which stores
// them in the ideals of gens, and each product is replaced by a small representative.
// The product is built from the last generator to the first, and the partial products
//      sufs[j] = a small representative of the class of &*[gens[i]^exps[i] : i in [j..#gens]]
// of the last call are stored in the last ideal of gens (attribute small_prod_cache): a call whose
// exponents agree with the previous ones from some position on starts from there. Magma enumerates
// the elements of an abelian group with the first coordinate changing fastest, so computing a
// representative for every class of a Picard group (as ICM does) costs about one product per class.
// The output might be a stored ideal (for example gens[i], if the only non-zero exponent is exps[i]
// = 1).
    assert #gens eq #exps;
    n:=#gens;
    if forall{ e : e in exps | e eq 0 } then
        return OneIdeal(S);
    end if;
    sufs:=[ OneIdeal(S) : j in [1..n+1] ]; // sufs[n+1] is the empty product
    j0:=n+1; // sufs[j] is known for j ge j0
    if assigned gens[n]`small_prod_cache then
        cgens,cexps,csufs:=Explode(gens[n]`small_prod_cache);
        if cgens eq gens then
            while j0 gt 1 and exps[j0-1] eq cexps[j0-1] do
                j0-:=1;
            end while;
            for j in [j0..n] do
                sufs[j]:=csufs[j];
            end for;
        end if;
    end if;
    for j in [j0-1..1 by -1] do
        e:=exps[j];
        if e eq 0 then
            sufs[j]:=sufs[j+1];
        else
            B:=small_power(gens[j],e);
            if forall{ i : i in [j+1..n] | exps[i] eq 0 } then
                sufs[j]:=B;
            else
                sufs[j]:=SmallRepresentative(sufs[j+1]*B);
            end if;
        end if;
    end for;
    last:=gens[n];
    last`small_prod_cache:=<gens,exps,sufs>;
    res:=sufs[1];
    res`IsInvertible:=true;
    return res;
end function;

//------------
// The maximal order
//------------

principal_by_short_element:=function(I)
// The cheap test, which often succeeds: whether I has a single generator, or is generated by
// ShortElement(I). If so, returns true and the generator; otherwise false, and I may still be
// principal. Nothing is stored: the intrinsic IsPrincipal stores the generator in I`Generators.
    if #Generators(I) eq 1 then
        return true,Generators(I)[1];
    end if;
    elt:=ShortElement(I);
    if elt*Order(I) eq I then
        return true,elt;
    end if;
    return false,_;
end function;

is_principal_maximal_components:=function(II, GRH)
// Given a fractional ideal II of the maximal order, returns whether it is principal and, if so, a
// generator (not stored, see IsPrincipal). The test is done in each component of the algebra, for a
// small representative of the class of II.
    assert IsMaximal(Order(II)); // this function should be called only for ideals of the maximal
                                 // order
    I,a:=SmallRepresentative(II); //a*II=I
    vprintf AlgEtQPicardGroup, 3:"is_principal_maximal_components:\n
                                    ZBasis(II) = %o\n,ZBasis(I) = %o\n",
                                    PrintSeqAlgEtQElt(ZBasis(II)),PrintSeqAlgEtQElt(ZBasis(I));
    test,I_asProd:=IsProductOfIdeals(I);
    assert test; // this function should be called only for ideals of the maximal order, hence I is
                 // a product
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
        // The class group is computed first: to prevent a bug of the in-built function IsPrincipal
        // (which might have been corrected by now...), and to use the GRH bound, since IsPrincipal
        // does not accept varargs.
        OL,oL:=class_group_internal(OrdIL, proof);
        testL,genL:=IsPrincipal(IL);
        assert2 (Zero(OL) eq (IL@@oL)) eq testL;
        if not testL then
            return false,_;
        end if;
        gen:=gen+embs[i](nf[i] ! genL);
    end for;
    assert2 gen*S eq I;

    gen:=gen/a;
    vprintf AlgEtQPicardGroup, 3:"is_principal_maximal_components:\n
                            [gen,a] = %o\n",PrintSeqAlgEtQElt([gen,a]);
    return true,gen;
end function;

is_principal_maximal:=function(II, GRH)
// Given a fractional ideal II of the maximal order, returns whether it is principal and, if so, a
// generator: first the cheap test, then the test in the components.
    test,gen:=principal_by_short_element(II);
    if test then
        return true,gen;
    end if;
    return is_principal_maximal_components(II,GRH);
end function;

picard_group_maximal:=function(O, GRH)
// Given the maximal order O, returns the Picard group of O, which is the direct sum of the class
// groups of the components, and the map from it to small representatives of the ideal classes, with
// preimages.
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
        GL,gL:=class_group_internal(OL, proof);
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
    end if;
    // the components of the unit ideal span O; each one alone does not span a fractional ideal, so
    // they are put together in one ideal (a single zero-divisor is not accepted by Ideal: for a
    // component Q)
    assert2 Ideal(O,&cat[ [ T[4](y) : y in Basis(T[2](Zero(T[1])),T[3]) ] : T in groups_maps_fields_maps ])
                                                                                        eq OneIdeal(O);
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
        // a small representative, which is invertible, hence 2-generated
        gen_O:=two_generated_invertible(SmallRepresentative(Ideal(O,gens_inA)),O);
        Append(~geninO,gen_O);
    end for;
    assert #geninO eq #Generators(G);
    rep_idinA:= function(x)
        vprint AlgEtQPicardGroup, 2: "picard_group_maximal: rep_idinA\n";
        coeff:=Eltseq(x);
        return small_prod_of_powers(O, geninO, coeff);
    end function;
    inverse_map:=function(id)
        vprint AlgEtQPicardGroup, 2: "picard_group_maximal: inverse_map\n";
        assert2 IsInvertible(id);
        // same ideal class, but integral and with a small Z-basis
        id:=SmallRepresentative(id);
        test,id_asprod:=IsProductOfIdeals(id);
        assert test;
        return &+[g[i](id_asprod[i]@@groups_maps_fields_maps[i,2]) : i in [1..#id_asprod]];
    end function;
    Codomain:=Parent(rep_idinA(Zero(G)));
    mapGtoO:=map<G -> Codomain | rep:-> rep_idinA(rep) , y:->inverse_map(y) >;
    assert2 forall{a : a in Generators(G)| (mapGtoO(a))@@mapGtoO eq a};
    vprintf AlgEtQPicardGroup, 2:"picard_group_maximal:\n
                          geninO = %o\n",[PrintSeqAlgEtQElt(ZBasis(I)) : I in geninO];
    return G,mapGtoO;
end function;

unit_group_maximal:=function(O, GRH)
// Given the maximal order O, returns its unit group, the direct sum of the unit groups of the
// components, and the map to the algebra, with preimages.
    assert IsMaximal(O); //this function should be used only for maximal orders
    test,O_asProd:=IsProductOfOrders(O);
    assert test; //the order must be a product
    A:=Algebra(O);
    nf,embs:=Components(A);
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
        vprint AlgEtQPicardGroup, 2: "unit_group_maximal: rep_inA\n";
        coeff:=Eltseq(rep);
        return &*[gensinA[i]^coeff[i] : i in [1..#coeff]];
    end function;

    disc_log:=function(x)
        vprint AlgEtQPicardGroup, 2: "unit_group_maximal: disc_log\n";
        comp_x:=Components(A ! x);
        x_in_Udp:=&*[ udp[i](comp_x[i]@@u_asProd[i]) : i in [1..#comp_x] ];
        return x_in_Udp;
    end function;

    maptoA:=map<Udp -> A | rep :-> rep_inA(rep) , y :-> disc_log(y) >;
    vprintf AlgEtQPicardGroup, 2:"unit_group_maximal:\n
                                    gensinA = %o\n",PrintSeqAlgEtQElt(gensinA);
    return Udp,maptoA;
end function;

//------------
// The exact sequence of a non-maximal order
//------------

// For a non-maximal order S with conductor F and maximal order O, with FO = F*O, the data of the
// exact sequence
//      1 -> S^* -> O^* --h--> B = (O/FO)^*/(S/F)^* -> Pic(S) -> Pic(O) -> 1,
// computed once by conductor_sequence and shared by UnitGroup and IsPrincipal (and by
// TotRealTotPos.m). PicardGroup reads R, (S/F)^* and UO from the same stored results, but not the
// record: see picard_group_non_maximal.
//  - O, F, FO; R, r: (O/FO)^* = ResidueRingUnits(O,FO); Sgens:
//    ResidueRingUnitsSubgroupGenerators(F), elements of S whose classes generate (S/F)^*;
//  - UO, uO: UnitGroup(O); B, b: the quotient of R by the classes of Sgens; h : UO -> B; K: the
//    kernel of h, which is S^* (UnitGroup(S) returns it);
//  - for reduce_mod_units_of_S: free, the indices of the free generators of UO; Kfree, the free
//    generators of K, as elements of UO; Klat (only if free is not empty), the lattice of their
//    exponents on the free generators of UO, with an LLL-reduced basis (see closest_vector_data);
//  - GRH: whether GRH was assumed.
conductor_sequence_fmt:=recformat< O, F, FO, R, r, Sgens, UO, uO, B, b, h, K, free, Kfree, Klat, GRH >;

conductor_sequence:=function(S, GRH)
// Returns the record conductor_sequence_fmt of the non-maximal order S, stored in S.
    if not assigned S`conductor_sequence or not valid_for(S`conductor_sequence`GRH,GRH) then
        // UO before the residue ring units, as UnitGroup(S) always did: both can make random
        // choices, so this order fixes the generators returned by UnitGroup(S)
        A:=Algebra(S);
        O:=MaximalOrder(A);
        UO,uO:=UnitGroup(O : GRH:=GRH );
        F:=Conductor(S);
        FO:=O!!F;
        R,r:=ResidueRingUnits(O,FO);
        Sgens:=ResidueRingUnitsSubgroupGenerators(F);
        B,b:=quo<R | [ a@@r : a in Sgens ]>;
        h:=hom<UO -> B | [ b(uO(UO.i)@@r) : i in [1..Ngens(UO)] ]>;
        K:=Kernel(h);
        data:=rec< conductor_sequence_fmt | O:=O, F:=F, FO:=FO, R:=R, r:=r, Sgens:=Sgens, UO:=UO, uO:=uO, B:=B,
                                            b:=b, h:=h, K:=K, GRH:=GRH >;
        data`free:=[ i : i in [1..Ngens(UO)] | Order(UO.i) eq 0 ];
        data`Kfree:=[ UO!K.i : i in [1..Ngens(K)] | Order(K.i) eq 0 ];
        assert #data`Kfree eq #data`free; // K has finite index in UO
        if #data`free gt 0 then
            data`Klat:=closest_vector_data(Matrix(Integers(),[ [ Eltseq(k)[i] : i in data`free ] :
                                                                                    k in data`Kfree ]));
        end if;
        S`conductor_sequence:=data;
    end if;
    return S`conductor_sequence;
end function;

reduce_mod_units_of_S:=function(u, data)
// Given u in UO and the record data = conductor_sequence(S), returns u-k with k in K = Kernel(h)
// (the units of S) such that the exponents of u-k on the free generators of UO are small: k is the
// element of K whose exponents are a closest vector to the ones of u.
    if #data`free eq 0 then
        return u;
    end if;
    v:=Vector([ Eltseq(u)[i] : i in data`free ]);
    c:=closest_vector_coefficients(data`Klat,v);
    k:=&+[ c[i]*data`Kfree[i] : i in [1..#c] ];
    assert2 data`h(k) eq Zero(data`B);
    return u-k;
end function;

is_principal_non_maximal:=function(I1, GRH)
// Given an invertible fractional ideal I1 of a non-maximal order S, returns whether it is principal
// and, if so, a generator (not stored, see IsPrincipal). The cheap test
// (principal_by_short_element) is not repeated.
    S:=Order(I1);
    A:=Algebra(S);
    O:=MaximalOrder(A);
    F:=Conductor(S);
    FO:=O!!F;
    cop,I:=CoprimeRepresentative(I1,F);
    IO:=O!!I;
    is_princ_IO,gen_IO:=is_principal_maximal(IO,GRH);
    if not is_princ_IO then
        return false,_;
    end if;
    IO`Generators:=[gen_IO]; // as IsPrincipal does
    R,r:=ResidueRingUnits(O,FO);
    if Order(R) eq 1 then
        // (O/FO)^* is trivial, so every generator of IO generates I
        assert2 gen_IO*S eq I;
        return true,gen_IO*cop^-1;
    end if;
    // I = cop*I1 is integral and coprime to F, so gen_IO, which lies in IO = I*O, is in O and is a
    // unit modulo FO. Every generator of I generates IO, so it is x = gen_IO/u for some u in O^*.
    // Conversely, such an x generates I iff its class in (O/FO)^* lies in (S/F)^*: at a prime P of
    // S not containing F we have S_P = O_P, so x*S_P = x*O_P = I_P; at a prime P with F c P we have
    // I_P = S_P, and x*S_P = S_P iff x is in S_P^*, iff x is in S_P (x is a unit of O_P, and S_P^*
    // = S_P meet O_P^* because O_P is integral over S_P), iff x mod F_P is in S_P/F_P. Writing B =
    // (O/FO)^*/(S/F)^* additively, the class of x in B is b(gen_IO) - h(u), where h : O^* -> B is
    // induced by the reduction modulo FO. Hence I is principal iff b(gen_IO) = h(u) for some u in
    // O^*, and then gen_IO/u is a generator of I.
    data:=conductor_sequence(S,GRH);
    is_princ,u:=HasPreimage(data`b(gen_IO@@r),data`h);
    if not is_princ then
        return false,_;
    end if;
    // u is determined modulo Kernel(h) = S^*. An arbitrary preimage can have huge exponents,
    // which makes the unit data`uO(u) (a product of powers of the fundamental units of O) huge.
    // The next call prevents this.
    u:=reduce_mod_units_of_S(u,data);
    gen_I:=gen_IO*(A ! data`uO(u))^-1;
    gen_I1:=gen_I*cop^-1;
    assert2 gen_I1*S eq I1;
    return true,gen_I1;
end function;

unit_group_non_maximal:=function(S, GRH)
// Given a non-maximal order S, returns its unit group, the kernel K of h : O^* -> B (see
// conductor_sequence), and the map to the algebra, with preimages.
    data:=conductor_sequence(S,GRH);
    UO:=data`UO; uO:=data`uO; B:=data`B; m:=data`h;
    P:=data`K;
    // The units of S are computed only when the map is evaluated, from their exponents on the
    // generators of O^*. When the index [O^*:S^*] is huge, these exponents are huge, already for
    // the generators of P, and such a unit cannot be written down as an element of A: for the
    // equation order of (x^5-x-1)(x^4+x^2+529)(x^3-2000) the generators of P have exponents up to
    // 10^25 (10^10 after an LLL reduction). So UnitGroup(S) computes only P, which is enough for
    // the inverse map and for the structure of S^*.
    map_P_to_S:=function(rep)
        vprint AlgEtQPicardGroup, 2: "UnitGroup: map_P_to_S\n";
        u:=UO!rep;
        // Magma computes x^n only for |n| < 2^30: a power of a non-torsion unit with such an
        // exponent would have billions of digits anyway.
        max_exp:=Max([ Abs(c) : c in Eltseq(u) ]);
        if max_exp ge 2^30 then
            error Sprintf("UnitGroup: this unit of the order is too large to be computed: its exponents on the generators of the unit group of the maximal order have up to %o digits. This happens when the index of the unit group of the order in the one of the maximal order is huge.", #Sprint(max_exp));
        end if;
        elt:=uO(u);
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
    p:=map<P -> Algebra(S) | x:->map_P_to_S(x), y:->map_S_to_P(y)  >;
    vprintf AlgEtQPicardGroup, 2:"UnitGroup:\n
                                exponents of the generators on the ones of O^* = %o\n",[ Eltseq(UO!P.i) : i in [1..Ngens(P)] ];
    return P,p;
end function;


//------------
// The Picard group of a non-maximal order
//------------

// The data of the Picard group of a non-maximal order S, computed by picard_group_non_maximal and
// used by the representatives (picard_representative) and the discrete logarithm
// (picard_discrete_log):
//  - S, its maximal order O, its conductor F, FO = F*O, and the map r of (O/FO)^* = R;
//  - GO, gO: PicardGroup(O), or the trivial group with no generators if Pic(O) is trivial (then gO
//    is not used); gens_GO_in_O: ideals of O, coprime to FO, in the classes of the generators of
//    GO; gens_GO_in_S: the invertible S-ideals gens_GO_in_O meet S, with two generators;
//  - H: the free abelian group on the generators of GO; D = R + H, with the injections mRD, mHD and
//    the projections mDR, mDH;
//  - P, mDP: Pic(S), the quotient of D by the relations: the classes of (S/F)^*, the images of O^*,
//    and the glue relations (the classes of the generators of GO, to the power of their orders, are
//    principal);
//  - generators_ideals: small representatives of the classes of the generators of P, with two
//    generators;
//  - GRH: whether GRH was assumed.
picard_fmt:=recformat< S, O, F, FO, r, GO, gO, gens_GO_in_S, gens_GO_in_O, H, D, mRD, mHD, mDR, mDH, P, mDP,
                       generators_ideals, GRH >;

picard_glue_generators:=function(S, O, FO, GO, gO)
// For every generator of GO = Pic(O): an ideal cI of O in its class, coprime to FO (see
// CoprimeRepresentative), and the invertible S-ideal cI meet S, with two generators. Returns the
// two sequences.
    gens_GO_in_S:=[]; //coprime with FO, in S and then meet S
    gens_GO_in_O:=[]; //coprime with FO, in O
    for i in [1..#Generators(GO)] do
        I:=gO(GO.i); //already made small in picard_group_maximal
        c,cI:=CoprimeRepresentative(I,FO);
        Append(~gens_GO_in_S,two_generated_invertible((S!!cI) meet S,S));
        Append(~gens_GO_in_O,cI); //used in the glue relations
    end for;
    return gens_GO_in_S,gens_GO_in_O;
end function;

picard_small_generators:=function(D, mDH, P, mDP)
// Returns preimages in D of the generators of P = D/Kernel(mDP), with small coefficients. The
// coefficients of a preimage on the free part H can be very large, which gives representatives with
// horrible Z-bases: each preimage g is moved by an element of the free part of the kernel close to
// it (see closest_vector_coefficients), the coefficients c on the generators of finite order N are
// replaced by c-N when this is smaller, and g itself is kept if it was smaller.
    gens_P_in_D:=[P.i@@mDP : i in [1..#Generators(P)]];
    ker:=Kernel(mDP);
    ker_gens_inD:=[D!ker.i:i in [1..Ngens(ker)]];
    ker_free_gens_inD:=[g:g in ker_gens_inD|Order(g) eq 0];
    if #ker_free_gens_inD eq 0 then
        return gens_P_in_D;
    end if;
    basis_mat:=Matrix([Eltseq(mDH(g)):g in ker_free_gens_inD]);
    lat:=closest_vector_data(basis_mat);
    vprintf AlgEtQPicardGroup, 2: "PicardGroup: basis_mat = \n%o\n",basis_mat;
    orders_Di:=[Order(D.i) : i in [1..Ngens(D)]];
    gens_P_in_D_small:=[];
    for g in gens_P_in_D do
        elt:=Vector(Eltseq(-mDH(g)));
        vprintf AlgEtQPicardGroup, 2: "PicardGroup: ClosestVector to %o...",elt;
        c:=closest_vector_coefficients(lat,elt);
        vprintf AlgEtQPicardGroup, 2: "done\n";
        g_small:=g+&+[c[i]*ker_free_gens_inD[i]:i in [1..#c]];
        // replace coeffs c of finite order N with N-c if better size
        c_small:=Eltseq(g_small);
        N_c_small:=[c_small[i]-orders_Di[i]:i in [1..#c_small]];
        c_better:=[orders_Di[i] ne 0 and (Abs(N_c_small[i]) lt Abs(c_small[i])) select N_c_small[i] else c_small[i]
                                                                                        : i in [1..#c_small]];
        g_small:=&+[c_better[i]*D.i:i in[1..Ngens(D)]];
        assert g-g_small in ker;
        vprintf AlgEtQPicardGroup, 2: "PicardGroup: g,g_small = %o\t%o\n",g,g_small;
        if &+[Abs(c):c in Eltseq(g)] lt &+[Abs(c):c in Eltseq(g_small)] then
            Append(~gens_P_in_D_small,g);
        else
            Append(~gens_P_in_D_small,g_small);
        end if;
    end for;
    return gens_P_in_D_small;
end function;

picard_representative:=function(data, rep)
// Given the record data of picard_fmt and an element rep of P, returns a small representative of
// the class rep.
    vprint AlgEtQPicardGroup, 3: "PicardGroup: picard_representative\n";
    return small_prod_of_powers(data`S, data`generators_ideals, Eltseq(rep));
end function;

compact_rep:=function(X, FO, r)
// Given an invertible ideal X of the maximal order O, coprime to FO, returns K = a*X, small,
// integral and coprime to FO, and the class of a in R = (O/FO)^* (the map r). Since X and K are
// coprime to FO, a is a unit at the primes above FO, so its class is defined. With MakeSmall, a is
// a short element of X^-1 also when X is already integral and coprime to FO, as the products in
// compact_power are.
    a,K:=CoprimeRepresentative(X,FO : MakeSmall:=true);
    return K,a@@r;
end function;

compact_power:=function(g, e, FO, r)
// Given an invertible ideal g of O coprime to FO and an integer e > 0, returns K and rho with g^e =
// mu^-1*K, K small, integral and coprime to FO, and rho the class of mu in R = (O/FO)^*. The power
// is never formed: repeated squaring, with every intermediate ideal replaced by compact_rep. The
// powers met along the way are stored in g`compact_powers, for the next calls (the discrete
// logarithms use the same powers again and again), together with FO.
    if not assigned g`compact_powers or not IsIdentical(g`compact_powers[1],FO) then
        g`compact_powers:=<FO,AssociativeArray()>;
    end if;
    is_def,Kr:=IsDefined(g`compact_powers[2],e);
    if is_def then
        return Explode(Kr);
    end if;
    if e eq 1 then
        K,rho:=compact_rep(g,FO,r);
    elif IsEven(e) then
        K1,r1:=$$(g,e div 2,FO,r);
        K,ra:=compact_rep(K1*K1,FO,r); // K = a*K1^2 = a*mu1^2*g^e
        rho:=2*r1+ra;
    else
        K1,r1:=$$(g,e-1,FO,r);
        K2,r2:=$$(g,1,FO,r);
        K,ra:=compact_rep(K1*K2,FO,r); // K = a*K1*K2 = a*mu1*mu2*g^e
        rho:=r1+r2+ra;
    end if;
    cache:=g`compact_powers[2];
    cache[e]:=<K,rho>;
    g`compact_powers:=<FO,cache>;
    return K,rho;
end function;

compact_power_product:=function(start, gens, exps, FO, r)
// Given invertible ideals start and gens[i] of O, coprime to FO, and integers exps[i] >= 0, returns
// K and rho with start*&*[gens[i]^exps[i]] = mu^-1*K, K small, integral and coprime to FO, and rho
// the class of mu in R = (O/FO)^*: see compact_power.
    K,rho:=compact_rep(start,FO,r);
    for i in [1..#gens] do
        if exps[i] ne 0 then
            Ki,ri:=compact_power(gens[i],exps[i],FO,r);
            K,ra:=compact_rep(K*Ki,FO,r);
            rho:=rho+ri+ra;
        end if;
    end for;
    return K,rho;
end function;

picard_discrete_log:=function(data, id)
// Given the record data of picard_fmt and an invertible ideal id of S, returns the element of P
// which is its class.
    vprint AlgEtQPicardGroup, 3: "PicardGroup: picard_discrete_log\n";
    O:=data`O; F:=data`F; FO:=data`FO; r:=data`r;
    // replace id with a small representative of the same class: otherwise every step below
    // (inverses, coprime representatives, principality tests) is very slow for ideals with large
    // Z-bases
    id:=SmallRepresentative(id);
    crep:=1/(CoprimeRepresentative(id^-1,F)); // here we check if id is invertible; (crep*id)^-1 is
                                              // coprime with F
    idO:=O!!(crep*id); //idO is coprime with FO
    GOrep:=idO@@data`gO;
    // prod = idO^-1 * &*[ gens_GO_in_O[i]^GOrep[i] ] is principal and coprime to FO (no negative
    // exponents): its generator elt gives the class of idO in Pic(S). The product is never formed:
    // prod = mu^-1*K with K small (see compact_power_product), and elt = gK/mu with gK a generator
    // of K, so the class of elt in R = (O/FO)^* is the one of gK minus rho, the class of mu.
    K,rho:=compact_power_product(idO^-1,data`gens_GO_in_O,Eltseq(data`H!Eltseq(GOrep)),FO,r);
    isprinc,gK:=is_principal_maximal(K,data`GRH);
    assert isprinc;
    Rrep:=gK@@r-rho;
    return data`mDP(-data`mRD(Rrep)+data`mHD(data`H!Eltseq(GOrep))); //[I]=-[xO meet S]+GOrep
end function;

picard_group_non_maximal:=function(S, GRH)
// Given a non-maximal order S, returns its Picard group, from the exact sequence (see
// conductor_sequence), and the map from it to small representatives of the ideal classes, with
// preimages.
    A:=Algebra(S);
    O:=MaximalOrder(A);
    GO,gO:=PicardGroup(O : GRH:=GRH); //stored in O
    F:=Conductor(S);
    FO:=O!!F;
    if #GO gt 1 then
        gens_GO_in_S,gens_GO_in_O:=picard_glue_generators(S,O,FO,GO,gO);
    else
        GO:=FreeAbelianGroup(0);
        gens_GO_in_S:=[];
        gens_GO_in_O:=[];
    end if;
    vprintf AlgEtQPicardGroup, 2:"PicardGroup:\n
                        gens_GO_in_S = %o\n
                        gens_GO_in_O = %o\n",
                        [PrintSeqAlgEtQElt(ZBasis(I)) : I in gens_GO_in_S],
                        [PrintSeqAlgEtQElt(ZBasis(I)) : I in gens_GO_in_O];

    // R, the generators of (S/F)^* and the units of O, from the stored results, in the order
    // PicardGroup always used: Magma can make random choices in these computations and in the
    // discrete logarithms in R, so this order fixes the generators of P below. For the same reason
    // the record conductor_sequence(S), which also computes the discrete logarithms of the units of
    // O in R, is not used here.
    R,r:=ResidueRingUnits(O,FO);
    Sgens:=ResidueRingUnitsSubgroupGenerators(F);
    UO,uO:=UnitGroup(O : GRH:=GRH );
    H:=FreeAbelianGroup(#Generators(GO));
    D, mRD, mHD, mDR, mDH := DirectSum(R,H);
    relDresidue:=[mRD(x@@r) : x in Sgens];
    relDunits:=[mRD(uO(x)@@r)  : x in Generators(UO)];
    // glue R and GO
    relDglue := [];
    assert #gens_GO_in_S eq #InvariantFactors(GO);
    for i in [1..#gens_GO_in_S] do
        n:=InvariantFactors(GO)[i];
        // g^n is principal, for g = gens_GO_in_O[i]: g^n = mu^-1*K with K small (see compact_power,
        // the power is never formed), and a generator of g^n is gK/mu with gK a generator of K,
        // whose class in R is the one of gK minus rho, the class of mu.
        K,rho:=compact_power(gens_GO_in_O[i],n,FO,r);
        is_princ, gK := is_principal_maximal(K,GRH);
        assert is_princ;
        Append(~relDglue,-mRD(gK@@r-rho)+mHD(H.i*n));
    end for;
    P, mDP := quo<D|relDresidue cat relDunits cat relDglue>;
    vprintf AlgEtQPicardGroup, 2: "PicardGroup: D = %o\n",D;
    vprintf AlgEtQPicardGroup, 2: "PicardGroup: P = %o\n",P;
    vprintf AlgEtQPicardGroup, 2: "PicardGroup: ker(mDP) = %o\n",Kernel(mDP);
    gens_P_in_D:=picard_small_generators(D,mDH,P,mDP);
    if #P eq 1 then
        return P,map<P->Parent(OneIdeal(S)) | rep:->OneIdeal(S), id:->Zero(P)>;
    end if;

    data:=rec< picard_fmt | S:=S, O:=O, F:=F, FO:=FO, r:=r, GO:=GO, gO:=gO, gens_GO_in_S:=gens_GO_in_S,
                            gens_GO_in_O:=gens_GO_in_O, H:=H, D:=D, mRD:=mRD, mHD:=mHD, mDR:=mDR, mDH:=mDH,
                            P:=P, mDP:=mDP, GRH:=GRH >;
    // the generators of P as small invertible ideals: [x O meet S] times the power product of the
    // gens_GO_in_S
    generators_ideals:=[];
    for igen->gen in gens_P_in_D do
        vprintf AlgEtQPicardGroup, 2: "PicardGroup: gens_P_in_D[%o] = %o...\n",igen,gen;
        id1:=(S!!( O*(r(mDR(gen))) )) meet S;
        // only the ideal class matters
        id2:=small_prod_of_powers(S, gens_GO_in_S, Eltseq(mDH(gen)));
        // a small representative, which is invertible, hence 2-generated
        Append(~generators_ideals,two_generated_invertible(SmallRepresentative(id1*id2),S));
    end for;
    data`generators_ideals:=generators_ideals;
    cod:=Parent(picard_representative(data,Zero(P)));
    p:=map<P -> cod | rep:->picard_representative(data,rep), id:->picard_discrete_log(data,id) >;
    vprintf AlgEtQPicardGroup, 2:"PicardGroup: computation is concluded\n
                        generators_ideals = %o\n",
                        [PrintSeqAlgEtQElt(ZBasis(I)) : I in generators_ideals];
    return P,p;
end function;

//------------
// Intrinsics
//------------

/// Returns whether the given fractional ideal is a principal ideal (over its order of definition), and, if so, also a generator. The parameter `GRH` (default false) can be set to true to assume GRH.
intrinsic IsPrincipal(I1::AlgEtQIdl : GRH:=false )->BoolElt, AlgEtQElt
{Returns whether the given fractional ideal is a principal ideal (over its order of definition), and, if so, also a generator. The parameter GRH (default false) can be set to true to assume GRH.}
    if not IsInvertible(I1) then
        return false,_;
    end if;
    test,gen:=principal_by_short_element(I1);
    if not test then
        if IsMaximal(Order(I1)) then
            test,gen:=is_principal_maximal_components(I1,GRH);
        else
            test,gen:=is_principal_non_maximal(I1,GRH);
        end if;
    end if;
    if not test then
        return false,_;
    end if;
    // the generator replaces the generators of I1
    I1`Generators:=[gen];
    return true,gen;
end intrinsic;

/// Returns the Picard group $\mathrm{Pic}(S)$ of the order $S$, which is not required to be maximal, and a map from the group to a set of representatives of the ideal classes. The parameter `GRH` (default false) can be set to true to assume GRH. The output is stored in $S$: a result computed without GRH is returned also when GRH is assumed, while a result computed under GRH is computed again when it is asked without GRH.
intrinsic PicardGroup( S::AlgEtQOrd : GRH:=false ) -> GrpAb, Map
{Returns the Picard group of the order S, which is not required to be maximal, and a map from the group to a set of representatives of the ideal classes. The parameter GRH (default false) can be set to true to assume GRH. The output is stored in S: a result computed without GRH is returned also when GRH is assumed, while a result computed under GRH is computed again when it is asked without GRH.}
    if not assigned S`PicardGroup or not valid_for(S`PicardGroup[3],GRH) then
        if IsMaximal(S) then
            P,p:=picard_group_maximal(S,GRH);
        else
            P,p:=picard_group_non_maximal(S,GRH);
        end if;
        S`PicardGroup:=<P,p,GRH>;
    end if;
    return S`PicardGroup[1],S`PicardGroup[2];
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

/// Returns the unit group of the given order. The parameter `GRH` (default false) can be set to true to assume GRH. The output is stored in $S$: a result computed without GRH is returned also when GRH is assumed, while a result computed under GRH is computed again when it is asked without GRH.
intrinsic UnitGroup(S::AlgEtQOrd : GRH:=false ) -> GrpAb, Map
{Returns the unit group of the given order. The parameter GRH (default false) can be set to true to assume GRH. The output is stored in S: a result computed without GRH is returned also when GRH is assumed, while a result computed under GRH is computed again when it is asked without GRH.}
    if not assigned S`UnitGroup or not valid_for(S`UnitGroup[3],GRH) then
        if IsMaximal(S) then
            U,u:=unit_group_maximal(S,GRH);
        else
            U,u:=unit_group_non_maximal(S,GRH);
        end if;
        S`UnitGroup:=<U,u,GRH>;
    end if;
    return S`UnitGroup[1],S`UnitGroup[2];
end intrinsic;

/// Given two fractional ideals $I$ and $J$ over the same order, checks whether $I=x\cdot J$, for some unit $x$. If so, $x$ is also returned. The parameter `GRH` (default false) can be set to true determines if the test is done under GRH.
intrinsic IsIsomorphic(I::AlgEtQIdl, J::AlgEtQIdl : GRH:=false ) -> BoolElt, AlgEtQElt
{Given two fractional ideals I and J over the same order, checks whether I=x*J, for some unit x. If so, x is also returned. The parameter GRH (default false) can be set to true determines if the test is done under GRH.}
    test,CIJ:=IsWeakEquivalent(I,J); // if so I=(I:J)*J and (I:J) is invertible in its
                                     // MultiplicatorRing
    if test then
        test2,x:= IsPrincipal(CIJ : GRH:=GRH );
        if test2 then
            assert2 I eq x*J;
            return test2,x;
        else
            return false, _ ;
        end if;
    else
        return false , _ ;
    end if;
end intrinsic;


/* TESTS

    printf "### Testing PicardGroup and UnitGroup:";
    //AttachSpec("../spec");
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

    // IsPrincipal and UnitGroup for a non-maximal order: E = Z[2000^(1/3)], #Pic(E) = 24, [O^*:E^*]
    // = 40 (Kluners-Pauli). The principal ideals p(g)*p(-g) with their generators, and p(g) for g
    // ne 0.
    _<x>:=PolynomialRing(Integers());
    A:=EtaleAlgebra(x^3-2000);
    E:=EquationOrder(A);
    P,p:=PicardGroup(E : GRH:=true);
    assert #P eq 24;
    for i in [1..10] do
        g:=Random(P);
        J:=p(g)*p(-g);
        test,gen:=IsPrincipal(J : GRH:=true);
        assert test and gen*E eq J;
        assert IsPrincipal(p(g) : GRH:=true) eq (g eq Zero(P));
    end for;
    U,u:=UnitGroup(E : GRH:=true);
    UO,uO:=UnitGroup(MaximalOrder(A) : GRH:=true);
    assert forall{ i : i in [1..Ngens(U)] | u(U.i) in E and u(U.i)^-1 in E };
    assert #quo<UO | [ u(U.i)@@uO : i in [1..Ngens(U)] ]> eq 40;
    printf ".";

    // UnitGroup of an order whose unit group has huge index in the one of the maximal order:
    // [O^*:E^*] has 30 digits. UnitGroup(E) computes the group, and a unit only when the map is
    // evaluated: the torsion units and the discrete logarithm work, while the generators of
    // infinite order are too large to be written down.
    K:=NumberField(x^4+x^2+529);
    M:=NumberField(x^3-2000);
    A:=EtaleAlgebra([K,K,M]);
    E:=EquationOrder(A);
    U,u:=UnitGroup(E : GRH:=true);
    UO,uO:=UnitGroup(MaximalOrder(A) : GRH:=true);
    assert TorsionFreeRank(U) eq TorsionFreeRank(UO) and TorsionFreeRank(U) eq 3;
    assert #Sprint(Index(UO,U)) eq 30;
    for i in [ i : i in [1..Ngens(U)] | Order(U.i) ne 0 ] do
        y:=u(U.i);
        assert y in E and y^-1 in E and y@@u eq U.i;
    end for;
    assert (-One(A))@@u in U;
    i:=[ i : i in [1..Ngens(U)] | Order(U.i) eq 0 ][1];
    try
        _:=u(U.i);
        too_large:=false;
    catch err
        too_large:=Position(err`Object,"too large to be computed") gt 0;
    end try;
    assert too_large;
    printf ".";

    // A representative for every class, as ICM computes them: the powers of the generators and the
    // partial products are stored and reused from one call to the next. Then some classes in a
    // random order.
    A:=EtaleAlgebra((x^2+5)*(x^2+20));
    E:=EquationOrder(A);
    P,p:=PicardGroup(E : GRH:=true);
    assert #P eq 320;
    elts:=[ g : g in P ];
    reps:=[ p(g) : g in elts ];
    assert forall{ I : I in reps | Order(I) eq E and I subset OneIdeal(E) };
    assert forall{ i : i in [1..#elts by 8] | reps[i]@@p eq elts[i] };
    for i in [1..20] do
        g:=Random(P);
        assert p(g)@@p eq g;
    end for;
    printf ".";

    //AttachSpec("../spec");
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

    //AttachSpec("../spec");
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

    //AttachSpec("../spec");
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
    // with GRH, as for P: a group computed without GRH would be a different one (see PicardGroup)
    C,c:=PicardGroup(O : GRH:=true);
    ext:=ExtensionHomPicardGroups(E,O : GRH:=true);
    for i in [1..100] do
        printf ".";
        x:=Random(P); y:=Random(P);
        assert ext(x+y) eq ext(x)+ext(y);
    end for;
    // When IsPrincipal returns true, the generator is stored as the only generator of the ideal, in
    // all the cases: the cheap test, maximal orders, and non-maximal orders with (O/FO)^* trivial
    // or not (Z+2O in Q(sqrt(-7)), where 2 splits, has (O/2O)^* trivial).
    _<z>:=PolynomialRing(Integers());
    K7:=EtaleAlgebra(z^2-z+2);
    O7:=MaximalOrder(K7);
    S7:=Order([ One(K7) ] cat [ 2*b : b in ZBasis(O7) ]);
    assert #ResidueRingUnits(O7,O7!!Conductor(S7)) eq 1;
    for S in [ S7, O7, EquationOrder(EtaleAlgebra(z^3-2000)), MaximalOrder(EtaleAlgebra(z^3-2000)) ] do
        s:=1;
        for k in [1..8] do
            gs:=[];
            for j in [1..3] do
                cs:=[];
                for i in [1..AbsoluteDimension(Algebra(S))] do
                    s:=(48271*s) mod (2^31-1);
                    Append(~cs,(s mod 7)-3);
                end for;
                Append(~gs,DotProduct(cs,ZBasis(S)));
            end for;
            if IsZeroDivisor(gs[1]) then
                continue;
            end if;
            I:=Ideal(S,gs);
            test,g:=IsPrincipal(I : GRH:=true);
            if test then
                assert Generators(I) eq [g] and g*S eq I;
            end if;
        end for;
    end for;
    printf ".";

    // An algebra with a component Q and a non-trivial class group: PicardGroup of the maximal order
    // failed in an assert2, which made an ideal of O from the single basis element of the component
    // Q, a zero-divisor.
    _<z>:=PolynomialRing(Integers()); // x and y are elements of P above
    A:=EtaleAlgebra((z-1)*(z^2+5));
    P,p:=PicardGroup(MaximalOrder(A) : GRH:=true);
    assert Invariants(P) eq [2] and forall{ g : g in P | p(g)@@p eq g };
    P,p:=PicardGroup(EquationOrder(A) : GRH:=true);
    assert forall{ g : g in P | p(g)@@p eq g };
    printf ".";

    // The stored results of PicardGroup and UnitGroup remember GRH: a result computed under GRH is
    // computed again when it is asked without GRH, and a result computed without GRH is returned
    // also under GRH.
    A:=EtaleAlgebra(z^3-2000);
    E:=EquationOrder(A);
    P1:=PicardGroup(E : GRH:=true);
    P2:=PicardGroup(E);
    P3:=PicardGroup(E : GRH:=true);
    assert not IsIdentical(P1,P2) and IsIdentical(P2,P3) and Invariants(P1) eq Invariants(P2);
    U1:=UnitGroup(E : GRH:=true);
    U2:=UnitGroup(E);
    U3:=UnitGroup(E : GRH:=true);
    assert not IsIdentical(U1,U2) and IsIdentical(U2,U3) and Invariants(U1) eq Invariants(U2);
    printf ".";
    printf " all good!";
*/
