/* vim: set syntax=magma :*/

freeze;

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
// Copyright 2024, S. Marseglia
/////////////////////////////////////////////////////

declare verbose PowerBass, 1;

import "../AlgEtQ/Ord.m" : crQZ , crZQ , Columns , hnf , MatrixAtoQ , MatrixAtoZ , MatrixQtoA , meet_zbasis , inclusion_matrix;

/* 
    Let M be an S-module in the K-algebra V.
    Assume that V is a 'pure-power' of K, that is, V = K^s for some positive integer s.
    If S is a Bass order then there are elements e1,...,es in M and fractional S-ideals I1,...,Is such that
    M = I1*e1 + ... + Is*es (direct sum).
    
    Relabel the fractional ideals such that
    (I1:I1) c (I2:I2) c ... c (Is:Is) [ the inclusion don't need to be strict ]
    Set Si:=(Ii:Ii).
    Then the S-isomorphism class of M is uniquely determined by the the tuple of orders (S1,...,Ss) and 
    the isomorphism class of the 'Steinitz class' I1*...*Is, which is an Ss-ideal.
    This a theorem of Bass.

    In this package, we inlcude functionalities to determine
    - if an S-module M has a 'Bass decomposition'
    - the Bass decomposition
    - the Steinitz class
    - the module M' in 'standard form' S1*e1+...+Ss*es where Si are overorders of S for i=1,..,s-1 and Ss is a fractional ideal
    - all isomorphism classes
*/


// ------------------- //
// Type AlgEtQMod
// ------------------- //

declare attributes AlgEtQMod :
                             DirectSumRep, // a sequence of pairs <J,m> where J is a practional R ideal and m:Algebra(R)->UniverseAlgebra, such that the module equals J1*m1(One(R))+...+Js*ms(One(R)) where the sum is direct. Also (Ji:Ji) subseteq (Ji+1:Ji+1)  
                             SteinitzClass, // the ideal J1*J2*...*Js of the above decomposition
                             StandardDirectSumRep, // < seq,map >, where seq is a sequence of pairs <Si,ei> with ei an orthogonal idempotent of the UniverseAlgebra UA and Si is the multiplicator ring of Ji (from the direct sum decomposition) for i = 1,...,s-1 and Ss is the SteinitzClass. map:UA->UA sends the DirectSumRep into seq.
                             multiplicator_overorder; // the multiplicator ring of the module, that is, the biggest over-order S of R for which M is an S-module

is_pure_power_internal:=function(m)
// Given a m:K->V module returns whether V is isomorphic to K^r, for some r, and, 
// if so, also a sequence of elements of V which are the r distinct copies if 1_K in V
// and an isomorphism K^r->V, together with embeddings K->Kr and projections Kr->K.
    K:=Domain(m);
    V:=Codomain(m);
    dimK:=AbsoluteDimension(K);
    dimV:=AbsoluteDimension(V);
    if dimV mod dimK ne 0 then
        return false,_,_,_,_;
    end if;
    r:=dimV div dimK;
    Knf:=[ DefiningPolynomial(L) : L in Components(K) ];
    Vnf:=[ DefiningPolynomial(L) : L in Components(V) ];
    if not #Seqset(Knf) eq #Components(K) then
        error "The number fields defining K are not distinct";
    end if;
    if not Seqset(Knf) eq Seqset(Vnf) then
        error "There is no componentwise diagonal action of K on V";
    end if;
    ind:=[ Index(Knf,fV) : fV in Vnf]; // there is only one index in Knf for each fV since 
                                       // Knf is a product of distinct number fields
    if not forall{ i : i in [1..#Knf] | #[ j : j in ind | j eq i] eq r } then
        return false,_,_,_,_;
    end if;
    oneKs:=[];
    Kis_in_V:=[ [ j : j in [1..#Vnf] | Vnf[j] eq Knf[i] ] : i in [1..#Knf] ]; 
        // [ [all occurences of Ki ] : i in [1..#Components(K)] ];
    assert Seqset(&cat(Kis_in_V)) eq {1..#Vnf}; 
    for i in [1..r] do
        pos:=[ ii[i] : ii in Kis_in_V ];
        ith_oneK:=V!<j in pos select 1 else 0 : j in [1..#Vnf]>;
        Append(~oneKs,ith_oneK);
    end for;
    if not m(One(K)) eq &+(oneKs) then
        return false,_,_,_,_;
    else
        Kr,embs,projs:=DirectProduct([K:i in [1..r]]);
        imgs:=[ &+[m(projs[i](b))*oneKs[i] : i in [1..r]] : b in AbsoluteBasis(Kr) ];
        isom:=Hom(Kr,V, imgs : ComputeInverse:=true);
        assert forall{ b : b in AbsoluteBasis(Kr) | Inverse(isom)(isom(b)) eq b } and isom(One(Kr)) eq &+(oneKs);
        return true,oneKs,isom,embs,projs;
    end if;
end function;

// ------------------- //
// functions for equality testing.
// ------------------- //

has_same_Zspan_Q:=function(gensM,gensN)
// input: two sequence of Qvectors
// returns: whether the vectors generate the same Z-module
    M_mat:=Matrix(gensM);
    N_mat:=Matrix(gensN);
    d_M:=Denominator(M_mat);
    d_N:=Denominator(N_mat);
    if not d_M eq d_N then 
    	return false;
    else
	    M_matZ:=crQZ(d_M*M_mat); 
	    N_matZ:=crQZ(d_N*N_mat);
        return hnf(M_matZ) eq hnf(N_matZ);
    end if;
end function;

has_same_Zspan_AlgEtQ:=function(gensM,gensN)
// given generators in the Universal Algebra, we check if they generate the same Z-module
    abs:=AbsoluteBasis(Algebra(gensM[1]));
    return has_same_Zspan_Q(AbsoluteCoordinates(gensM,abs),AbsoluteCoordinates(gensN,abs));
end function;

// ------------------- //
// Multiplicator Ring
// ------------------- //

compute_multiplicator_overorder:=function(M)
// given an AlgEtQMod over the order R, it returns the biggest overorder S of R such that S*M=M
    if not assigned M`multiplicator_overorder then
        if assigned M`DirectSumRep then
            mult_rings:=[ MultiplicatorRing(J[1]) : J in M`DirectSumRep ];
            assert exists(mr){ T : T in mult_rings | forall{ S : S in mult_rings | T subset S }  }; //is the smallest of the multiplicator rings
            assert mr in mult_rings;
            M`multiplicator_overorder:=mr;
        else
            gens:=ZBasis(M);
            _,map:=UniverseAlgebra(M);
            R:=Order(M);
            oo:=OverOrders(R);
            oo1:=[ S : S in oo | has_same_Zspan_AlgEtQ( [ g*map(s) : s in ZBasis(S) , g in gens ] cat gens , gens ) ];
            mr:=Order(&cat[ZBasis(S) : S in oo1]);
            assert mr in oo;
            M`multiplicator_overorder:=mr;
        end if;
    end if;
    return M`multiplicator_overorder;
end function;

// ------------------- //
// DirectSumRep
// ------------------- //

bass_factor:=function(gens_input,Rg,mAgs,embs,projs)
// Input: 
// - gens_input is a sequence of vectors in mAgs:Ag->Ag^s whose Zspan is an Rg-module M, where Ag is the fraction field of Rg.
// - Rg : see above
// - mAgs: see above
// - embs,projs are the embeddings Ag->Ags and projections Ags->Ag
// Output: J,v,gens_M0 where 
// - J is a fractional Rg-ideal, 
// - v is an element of M c Ags
// - gens_M0 are the generators of a sub-Rg-module M0 of M 
// such that M = J*v + M0 (direct sum).

    Ag:=Algebra(Rg);
    Ags:=Algebra(gens_input[1]);
    n:=AbsoluteDimension(Ag);
    abs_Ags:=AbsoluteBasis(Ags);
    ns:=AbsoluteDimension(Ags);
    s:=ns div n;
    assert s eq #projs;

    // We search among OverOrder(Rg) for the biggest order such that TM = M.
    T:=Order(&cat [ZBasis(S):S in OverOrders(Rg)|has_same_Zspan_AlgEtQ([g*mAgs(s) :s in ZBasis(S), g in gens_input] cat gens_input,gens_input)]);

    // Let u be a unit of Ags. Write Ags = Ag*u + B where B is isomorphic to Ag^(s-1).
    // Let p_u be the projection along u defined by sending x*u+b (with x in Ag and b in B) to x.
    //
    // We need a generators u from gens such that th fractional Rg-ideal J:=p_u(M) satisfies (J:J)=T.
    // If no u in gens satisfies this property, we act with a matrix rU in GL_ns(Z) on gens, 
    // so that the obtained set also has Zspan M, until we find one.
    // An element u with the property above exists in M, by a Theorem of Bass.
    //
    // TODO It would be nice to find a deterministic method.
    k1:=1;
    tt:=false;
    gens_input_matrix:=MatrixAtoQ(gens_input); // matrix with Q entries
    //s0:=Rank(Matrix(gens_input));
    s0:=Rank(gens_input_matrix); // rank of M over Z
    rU:=IdentityMatrix(Rationals(),s0);
    repeat
        for k2 in [1..100] do
            gens:=MatrixQtoA(Ags,rU*gens_input_matrix);
            assert has_same_Zspan_AlgEtQ(gens,gens_input);
            for ui in [1..#gens] do
                u:=gens[ui];
                if not IsZeroDivisor(u) then
                    vprintf PowerBass,2: ".";
                    // We define the 'projection along u', as described above.
                    p_u:=[ projs[j](u) : j in [1..s] ]; // the projections of u onto the s Ag-components of Ags 
                    n_u:=&+[ (p_u[j])^2 : j in [1..s] ]; // the 'Ag-norm', squared, of u
                    p_u:=Hom(Ags,Ag,[(&+[ projs[j](g)*p_u[j]:j in [1..s]])/n_u:g in abs_Ags]);
                    J:=Ideal(Rg,[ p_u(g) : g in gens ]); // J=p_u(M)
                    S:=MultiplicatorRing(J);
                    assert T subset S;
                    if S eq T then
                        SwapElements(~gens,1,ui); // we put u in first position
                        tt:=true;
                        break k2;
                    end if;
                end if;
            end for;
            rU:=crZQ(RandomUnimodularMatrix(s0,10^k1));
        end for;
        k1+:=1;
    until tt; //tt=true if S eq T
    
    // We got our factor J. Need to find v and M0
    Jinv:=ColonIdeal(S,S!!J);
    gens_JinvM:=&cat[[ mAgs(j)*g :g in gens] : j in ZBasis(Jinv)]; // generators of J^-1*M (not a ZBasis)
    coeffsu_JinvM:=[ p_u(g) : g in gens_JinvM ]; // p_u(J^-1*M)
    assert Ideal(Rg,coeffsu_JinvM) eq Ideal(Rg,ZBasis(S));
    A:=MatrixAtoQ(coeffsu_JinvM);
    d_A:=Denominator(A);
    AZ:=crQZ(d_A*A); 
    assert #Rows(AZ) eq #gens_JinvM;
    vect:=Vector([ d_A ] cat [0 : i in [1..n-1] ]);
    v:=SumOfProducts(Eltseq(Solution(AZ,vect)),gens_JinvM);

    // Define M0=J^-1*N0, where N0 is the set of elements of J^-1*M such that the projection onto u is 0, 
    // that is, N0=ker(J^-1*M ->> Ag*u)    
    
    mat:=MatrixAtoQ([p_u(g)*n_u : g in gens_JinvM]);
    den:=Denominator(mat);
    mat:=crQZ(den*mat); // over Z
    s,ker:=Solution(mat,Vector([0 : i in [1..NumberOfColumns(mat)]]));
    assert IsZero(s); // Zero is always a solution, s needs to be Zero. 
                      // If that is not the case I need to modify the generators of N0 accordingly.
    // N0 is generated by ker over Z.    
    if Dimension(ker) eq 0 then 
        gens_M0:=[];
    else
        gens_N0:=MatrixQtoA(Ags,crZQ(Matrix(Basis(ker)))*MatrixAtoQ(gens_JinvM));
        // M0 = J*N0
        gens_M0:=[ mAgs(j)*m : m in gens_N0 , j in ZBasis(J) ];
        // reduce the number using hnf
        mat_gensZ_M0:=MatrixAtoQ(gens_M0);
        den:=Denominator(mat_gensZ_M0);
        gens_M0:=MatrixQtoA(Ags,(1/den)*crZQ(hnf(crQZ(den*mat_gensZ_M0))));
    end if;
    vprint PowerBass : J,v,gens_M0;
    // We make J smaller, and modify v accordingly.
    J,sm:=SmallRepresentative(J);
    v:=v/mAgs(sm);
    //TEST//
        gens_Jv:=[ mAgs(j)*v : j in ZBasis(J)];
        assert has_same_Zspan_AlgEtQ(gens_M0 cat gens_Jv,gens_input);
        assert #gens_M0 eq #gens - n;
        if not IsEmpty(gens_M0) then
            assert Rank(MatrixAtoQ(gens_M0)) eq Rank(MatrixAtoQ(gens_input))-n;
            assert Rank(Matrix([[p(g):p in projs]: g in gens_M0])) eq Rank(Matrix([[p(g):p in projs]: g in gens_input]))-1;
        else
            assert n eq Rank(MatrixAtoQ(gens_input));
            assert 1 eq Rank(Matrix([[p(g):p in projs]: g in gens_input]));
        end if;
    //END TEST//    
    return J,v,gens_M0;
end function;

intrinsic DirectSumRepresentation(M::AlgEtQMod)->SeqEnum[Tup]
{Given an R-module M in V=K^s, with R Bass, it returns a sequence of pairs <J,m> where J is a practional R ideal and m:Algebra(R)->UniverseAlgebra, such that M =  J1*m1(One(R))+...+Js*ms(One(R)) where the sum is direct. Also (Ji:Ji) subseteq (Ji+1:Ji+1).}
    if not assigned M`DirectSumRep then    
        UA,mR:=UniverseAlgebra(M);
        Rg:=Order(M);
        is_pure_power,ones,AgstoUA,embs,projs:=is_pure_power_internal(mR); //AGstoUA:Ags^s->UA isomorphism
        Ags:=Domain(AgstoUA);
        require is_pure_power : "The universe algebra is not of the form K^r.";
        Ag:=Algebra(Rg);
        n:=Dimension(Ag);
        s:=Dimension(Ags) div n;
        if IsBass(Rg) then
            T:=compute_multiplicator_overorder(M);
            UAtoAgs:=Inverse(AgstoUA);
            mAgs:=mR*UAtoAgs;
            gens:=[ UAtoAgs(g) : g in ZBasis(M) ]; // in Ags
            // We start working in Ags, instead of UA. We will go back at the end.
            gens_loop:=gens;
            Aloop:=Ags;
            embs_loop:=embs;
            projs_loop:=projs;
            mAg_Aloop:=mAgs;
            mAloop_Ags:=Hom(Aloop,Ags,AbsoluteBasis(Aloop)); //identity
            bass_decompI:=[];
            rank:=s;
            repeat
                vprintf PowerBass,2 : "\ngens_loop = %o",gens_loop;
                J,v0,gens_I0:=bass_factor(gens_loop,Rg,mAg_Aloop,embs_loop,projs_loop);
                rank-:=1; // rank=rk_Rg(<gens_I0>_Z)
                // v0 and gens_I0 are in Aloop
                // I need to move v0 to Ags
                v:=mAloop_Ags(v0);
                Append(~bass_decompI,<J,v>);
                assert #gens_I0 eq rank*n;
                // Assume that rk_Rg(<gens_loop>_Z) = r.
                // Then rk_Rg(<gens_I0>_Z) = r-1.
                // he next run of bass_factor should happen in Ag^(r-1), appropriately embedded is Ag^s,
                // where s=rk_Rg(M).
                if rank ge 1 then  
                    mI0:=Matrix([[p(g):p in projs_loop]:g in gens_I0]);
                    rkI0:=Rank(mI0);
                    rowsI0:=Rows(mI0);
                    repeat
                        m1_mat:=Matrix( [Random(rowsI0) : i in [1..rkI0]] );
                    until Rank(m1_mat) eq rkI0;
                    assert #Columns(m1_mat) eq #embs_loop;
                    m1:=[ &+[embs_loop[i](g[i]) : i in [1..#embs_loop]] : g in Rows(m1_mat) ]; // as elts of Aloop
                    Aold:=Aloop;
                    Aloop,embs_loop,projs_loop:=DirectProduct([Ag : i in [1..rank]]);
                    mAloop_Aold:=Hom(Aloop,Aold, 
                        [ &+[ mAg_Aloop(projs_loop[i](b))*m1[i] : i in [1..rank] ] : b in AbsoluteBasis(Aloop) ]);
                    assert forall{ i : i in [1..rank] | mAloop_Aold(ort[i]) eq m1[i] where ort:=OrthogonalIdempotents(Aloop)}; //needs to send the orthogonal idempotents to the elts in m1
                    mAloop_Ags:=mAloop_Aold*mAloop_Ags;
                    mAg_Aloop:=NaturalAction(Ag,Aloop);
                    gens_loop:=Rows(Solution(m1_mat,mI0)); //gens_loop*m1 eq mI0, as matrix
                    gens_loop:=[&+[embs_loop[i](g[i]):i in [1..rank]]:g in gens_loop]; // as elts of Aloop
                end if;
            until rank eq 0;
            assert #gens_I0 eq 0;

            // We need to push back all the v's from Ags to UA
            basisAg:=AbsoluteBasis(Ag);
            BC:= [ <bc[1],Hom(Ag,UA,[mR(basisAg[i])*AgstoUA(bc[2]): i in [1..Dimension(Ag)]]: ComputeInverse:=false)> : bc in bass_decompI ];
            // TESTS //////////////////////////////////
                BC_test:=[[ bc[2](z) : z in ZBasis(bc[1]) ] : bc in BC ];
                assert has_same_Zspan_AlgEtQ( ZBasis(M) , &cat(BC_test)); // the sum is ok
                assert forall{bc : bc in BC | Codomain(bc[2]) eq UA};
                // test using abelian groups. I will test also the intersections
                mat_gens:=MatrixAtoQ(ZBasis(M));
                F:=FreeAbelianGroup(Dimension(Ags));
                subs:=[];
                for gens_Jv in BC_test do
                    mat_gens_Jv:=MatrixAtoQ(gens_Jv);
                    den:=Denominator(VerticalJoin(mat_gens,mat_gens_Jv));
                    coords_Jv:=Solution(crQZ(den*mat_gens),crQZ(den*mat_gens_Jv));
                    gensF_Jv:=[F!Eltseq(r) : r in Rows(coords_Jv) ];
                    subF_Jv:=sub<F|gensF_Jv>;
                    Append(~subs,subF_Jv);
                end for;
                assert &+subs eq F; //the sum is ok
                assert forall{H : H in subs | not exists{G : G in subs | G ne H and not IsTrivial(G meet H)}}; 
                    //the intersections are trivial, that is the sum is direct
            // END TESTS //////////////////////////////////////
            M`DirectSumRep:=BC;
        elif n eq Dimension(UA) then //squarefree case
            M`DirectSumRep:=<Ideal(Rg,ZBasis(M)),mR>;
        else
            error "implemented only for the squarefree or power-of-bass case";
        end if;
    end if;
    return M`DirectSumRep;
end intrinsic;

// ------------------- //
// Compute all Modules up to isomorphism
// ------------------- //

chains_of_overorders:=procedure(r,~chain,~all_chains)
	if #chain eq r then 
		Append(~all_chains,chain); 
	else
		T:=chain[#chain];        
		for S in OverOrders(T) do
			chain1:=chain;
			Append(~chain1,S);
			$$(r,~chain1,~all_chains);
		end for;
	end if;
end procedure;

direct_sums_overorders:=function(R,r)
//given a (Bass) order R and a rank r, it returns all sequences S_1,..,S_r where R\subseteq S_i
	choo:=<>;
	oo:=OverOrders(R);
	for S in oo do
		ch:=<S>;
		chains_of_overorders(r,~ch,~choo);
	end for;
	return choo;
end function;

intrinsic IsomorphismClassesOverBassOrder(R::AlgEtQOrd, r::RngIntElt)->SeqEnum[AlgEtQMod]
{Given a Bass order R and a rank r, it returns representatives of all the isormorphism classes of torsion-free R modules of rank r.}
	require r gt 0 : "the second input must be a positive integer";  
    Kis:=Components(Algebra(R));
    V:=EtaleAlgebra([Kis[j] : i in [1..r] , j in [1..#Kis]]); // V=K^r
    map:=NaturalAction(Algebra(R),V);
	return IsomorphismClassesOverBassOrder(R,map);
end intrinsic;

intrinsic IsomorphismClassesOverBassOrder(R::AlgEtQOrd, map::Map )->SeqEnum[AlgEtQMod]
{Given a Bass order R and a map:Algebra(R)->UA, it returns representatives of all the isormorphism classes of torsion-free R modules inside R. We need UA to be isomorphic to a power of Algebra(R).}
    is_pure_power,ones:=is_pure_power_internal(map);
    require is_pure_power : "The universe algebra is not of the form K^r.";
    require IsBass(R) : "the first input must be a Bass order";
    UA:=Codomain(map);
    assert Type(UA) eq AlgEtQ;
    AR:=Algebra(R);
    assert Dimension(UA) mod Dimension(AR) eq 0;
    r:=Dimension(UA) div Dimension(AR);
    seqs:=<>;
    choo:=direct_sums_overorders(R,r);
    for ch in choo do
        P,p:=PicardGroup(ch[r]);
        pic:=[p(g) : g in P];
        ch_rem:=Prune(ch);
        for I in pic do
            Append(~seqs,Append(ch_rem,I));
        end for;
    end for;
    output:=[];
    for seq in seqs do
        assert #seq eq #ones;
        seqid:=[ Ideal(R,ZBasis(seq[i])) : i in [1..#seq] ];
        assert2 forall{ i : i in [2..#seqid] | MultiplicatorRing(seqid[i-1]) subset MultiplicatorRing(seqid[i])};
        // the conditions asserted in the previous condition is assumed by the DirectSumRep attribute
        M:=ModuleFromDirectSum(R,map,[<seqid[i],ones[i]> : i in [1..#seq]]);
        maps_ones:=[ Hom(AR,UA,[map(a)*ones[i] : a in AbsoluteBasis(AR)]) : i in [1..#ones] ];
        M`DirectSumRep:=[ <seqid[i],maps_ones[i]> : i in [1..#ones] ]; 
        Append(~output,M);
    end for;
    return output;
end intrinsic;

// ------------------- //
// SteinitzClass
// ------------------- //

intrinsic SteinitzClass(M::AlgEtQMod)->AlgEtQOrdIdl
{returns the product of the fractional ideals in the input}
    if not assigned M`SteinitzClass then
	    M`SteinitzClass:=&*[ D[1] : D in DirectSumRepresentation(M) ];
    end if;
	return M`SteinitzClass;
end intrinsic; 

// ------------------- //
// StandardDirectSumRep
// ------------------- //

std_Bass_map:=function(bc,UA,isom,embs,projs)
// given a sequence I1,...,Ir of ideals of K it returns the hom:UA->UA that 
// sends I1+...+Is to S1+...+Ss-1+I1*..*Is acting on column vectors
// with S1 Bass and S1 in S2 in S3 in ... in Ss-1 in Sr = (I1*...*Is:I1*...*Is)
    s:=#bc;
    AR:=Algebra(bc[1]);
    id:=IdentityMatrix(AR,s);
    M:=id;
    for i in [1..s-1] do
        I1:=bc[i];
        I2:=bc[i+1];
        S:=MultiplicatorRing(I1);
        c2I2,c2:=MakeIntegral(S!!I2);
        c2:=AR!c2;
        c1,c1I1:=CoprimeRepresentative(S!!I1,c2I2);
        a2:=ChineseRemainderTheorem(c1I1,c2I2,AR!1,AR!0);
        a1:=1-a2;
        Mi:=InsertBlock(id,Matrix(2,2,[c1,-c2,a2/c2,a1/c1]),i,i); //acts on column vectors
        M:=Mi*M;
        bc[i+1]:=I1*I2; //this was missing
    end for;
    MT:=Transpose(M);
    ARstoUA:=func< vec | isom(&+[embs[i](vec[i]) : i in [1..#vec]])>;
    UAtoARs:=func< x | [projs[i](Inverse(isom)(x)) : i in [1..#projs]]>;
    images:=[ ARstoUA(Rows(Matrix(1,s,UAtoARs(b))*MT)[1]) : b in Basis(UA) ];
    hom_id:=Hom(UA,UA,images : ComputeInverse:=true);
    return hom_id;
end function;

intrinsic StandardDirectSumRepresentation(M::AlgEtQMod)->SeqEnum[Tup],Map
{Returns seq,map, where seq is a sequence of pairs <Si,ei> with:
    - ei:Algebra(R)->UA that sends 1_R to the ith orthogonal idempotent of the UniverseAlgebra UA, and
    - Si is the multiplicator ring of Ji (from the DirectSumRepresentation) for i = 1,...,s-1 and Ss is the SteinitzClass.
      [ note : S1 in S2 in ... in Ss-1 in (Ss:Ss) ].
    - map:UA->UA is an R-linear morphism that sends the DirectSumRepresentation into seq.}
    if not assigned M`StandardDirectSumRep then
        UA,mR:=UniverseAlgebra(M);
        R:=Order(M);
        is_pure_power,oid,isom,embs,projs:=is_pure_power_internal(mR);
        require is_pure_power : "The universe algebra is not of the form K^r.";
        AR:=Algebra(R);
        oneAR:=One(AR);
        dr:=Dimension(AR);
        DRS:=DirectSumRepresentation(M);
        s:=#DRS;
        assert #oid eq s;
        Si_s:=[ i lt s select Ideal(R,ZBasis(MultiplicatorRing(DRS[i,1]))) else SteinitzClass(M) : i in [1..s] ]; 
        ei_s:=[ Hom(AR,UA,[mR(b)*oid[i] : b in Basis(AR) ] : ComputeInverse:=false) : i in [1..s] ];
        seq:=[ < Si_s[j],ei_s[j] > : j in [1..s ]];
        //now we build map:
        hom_id:=std_Bass_map([D[1] : D in DRS],UA,isom,embs,projs);
        hom_vs:=Inverse(Hom(UA,UA,&cat[ [mR(a)*D[2](oneAR):a in Basis(AR)] : D in DRS ] : ComputeInverse:=true));
        map:=Hom(UA,UA,[hom_id(hom_vs(b)):b in Basis(UA)]:ComputeInverse:=true); // not as a composition to make it cleaner.
        M`StandardDirectSumRep:=<seq,map>;
        // TEST 
        if GetAssertions() gt 1 then
            gens_in:=&cat[[ mR(z)*D[2](oneAR) : z in ZBasis(D[1])] : D in DRS];
            gens_std:=&cat[[ mR(z)*D[2](oneAR) : z in ZBasis(D[1])] : D in M`StandardDirectSumRep[1]];
            assert has_same_Zspan_AlgEtQ(ZBasis(M),gens_in);
            assert has_same_Zspan_AlgEtQ([ map(g) : g in gens_in],gens_std);
        end if;
        // end test
    end if;
	return M`StandardDirectSumRep[1],M`StandardDirectSumRep[2];
end intrinsic;

// ------------------- //
// Isomorphism testing for AlgEtQMod over Bass Orders
// ------------------- //

intrinsic IsIsomorphicOverBassOrder( M1::AlgEtQMod, M2::AlgEtQMod )->BoolElt,Map
{Returns wheter two AlgEtQMods M1 and M2 are isomorphic and if so it returns also a map from the common universe algebra that sends M1 into M2.}    
    UA,mR:=UniverseAlgebra(M1);
    R:=Order(M1);
    AR:=Algebra(R);
    require UniverseAlgebra(M2) eq UA : "the modules don't live in the same algebra";
    require R eq Order(M2) : "the modules are not defiend over the same order";
    is_pure_power,ones,isom,embs,projs:=is_pure_power_internal(mR);
    require is_pure_power : "The universe algebra is not of the form K^r.";
    Std1,mStd1:=StandardDirectSumRepresentation(M1);   
    Std2,mStd2:=StandardDirectSumRepresentation(M2);
    s:=#Std1;
    test:=&and[ Std1[j,1] eq Std2[j,1] : j in [1..s-1] ];
    if test then
        test,x:=IsIsomorphic(Std1[s,1],Std2[s,1]); // Std1[s]=x*Std2[s]
        if test then
            ARstoUA:=func< vec | isom(&+[embs[i](vec[i]) : i in [1..#vec]])>;
            UAtoARs:=func< x | [projs[i](Inverse(isom)(x)) : i in [1..#projs]]>;
            images:=[ARstoUA([ j lt s select UAtoARs(b)[j] else UAtoARs(b)[j]/x :j in [1..s]]):b in Basis(UA)];
            map_x:=Hom(UA,UA,images); 
            // TEST
            if GetAssertions() gt 1 then
                oneAR:=AR!1;
                gens1:=&cat[[ mR(z)*D[2](oneAR) : z in ZBasis(D[1])] : D in Std1 ];
                gens2:=&cat[[ mR(z)*D[2](oneAR) : z in ZBasis(D[1])] : D in Std2 ];
                assert has_same_Zspan_AlgEtQ([map_x(g) : g in gens1] , gens2);
            end if;
            // end test
            map:=Hom(UA,UA,[ Inverse(mStd2)(map_x(mStd1(g))) : g in Basis(UA) ]);
            // TEST
            assert has_same_Zspan_AlgEtQ([map(g) : g in ZBasis(M1)] , ZBasis(M2));
            // end test
            return true,map;
         else
            return false,_;
         end if;
    else
        return false,_;
    end if;
end intrinsic;       

/*
//TESTS

    //TODO move these tests to all_tests_AlgEtQMod.m

    AttachSpec("~/AlgEt/spec");
    AttachSpec("~/AlgEt/specMod");
    AttachSpec("~/AlgEt/specMtrx");
    SetDebugOnError(true);
    
    P<x>:=PolynomialRing(Integers());
    f:=x^2 + x + 2;
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    R:=Order([2*F]);

    O:=MaximalOrder(A);
    ff:=Conductor(R);
    V:=DirectProduct([A,A,A]);
    m:=NaturalAction(A,V);
    Vnf:=Components(V);
    Anfpoly:=[ DefiningPolynomial(L) : L in Components(A) ];
    Vnfpoly:=[ DefiningPolynomial(L) : L in Vnf ];
    MO:=Module(R,m,<1*MaximalOrder(Vnf[i]) : i in [1..#Vnf]>);
    ff:=O!!Conductor(R);
    test,ff_prod:=IsProductOfIdeals(ff);
    assert test;
    ind:=[ #Vnfpoly+1 - Index(Reverse(Vnfpoly),fA)  : fA in Anfpoly ]; // last occurence of each 
                                                                       // number field from the dec of K 
                                                                       // in the decomposition of V
    Mff:=Module(R,m,<ff_prod[Index(Anfpoly,Vnfpoly[i])] : i in [1..#Vnf]>);
    candidates:=MaximalIntermediateModules(MO,Mff);
    for M in candidates do
        _:=DirectSumRepresentation(M);
        _:=SteinitzClass(M);
    end for;

    ///////////////////

    AttachSpec("~/AlgEt/spec");
    AttachSpec("~/AlgEt/specMod");
    AttachSpec("~/AlgEt/specMtrx");
    SetDebugOnError(true);

    // takes a while
    P<x>:=PolynomialRing(Integers());
    f:=x^4 + 3*x^3 + 8*x^2 + 39*x + 169;
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    q:=Integers() ! (Coefficients(f)[1]^(2/Degree(f)));
    R:=Order([F,q/F]);
    ver:=[62,97,144,206,286];
    res:=[ #IsomorphismClassesOverBassOrder(R,i) : i in [1..5] ];
    assert res eq ver;
    
    ///////////////////
    
    AttachSpec("~/AlgEt/spec");
    AttachSpec("~/AlgEt/specMod");
    AttachSpec("~/AlgEt/specMtrx");
    SetDebugOnError(true);

    P<x>:=PolynomialRing(Integers());
    f:=x^6-x^5+2*x^4-2*x^3+4*x^2-4*x+8;
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    R:=Order([F,ComplexConjugate(F)]);
    iso_cl:=IsomorphismClassesOverBassOrder(R,3);
    assert #iso_cl eq 6;

    M:=iso_cl[1];
    
    _:=[ StandardDirectSumRepresentation(M) : M in iso_cl ];

    t,s:=IsIsomorphicOverBassOrder(iso_cl[1],iso_cl[1]);
    assert t;
    t,s:=IsIsomorphicOverBassOrder(iso_cl[1],iso_cl[2]);
    assert not t;

    for M1,M2 in iso_cl do
        t,s:=IsIsomorphicOverBassOrder(M1,M2);
        assert2 IsIsomorphic(M1,M2); // generic test
    end for;

*/
