/* vim: set syntax=magma :*/

freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

// declare verbose ?????, 1;

//------------
// 
//------------

import "../LowCohenMacaulayType.m" : wkicm_bar_CM_type2;

intrinsic IntermediateIdealsVSWithTrivialExtensionAndPrescribedMultiplicatorRing(I::AlgEtQIdl,P::AlgEtQIdl, O::AlgEtQOrd)->SetIndx[AlgEtQIdl]
{Given a fractional S-ideal I, a prime ideal P of S, and an overorder O of S such that O subset (I:I), it returns all the fractional S-ideal J such that PI subset J subset I such that (J:J) eq S and JO=I.}
    require IsPrime(P) : "The ideal P must be prime";
    S:=Order(P);
    require S eq Order(I) : "The ideals must be over the same order";
    require S subset O : "O is not an overorder of Order(I)";
    IO:=O!!I;
    require O subset MultiplicatorRing(I) : "I is not an O-ideal";
    PI:=P*I;
    zbPI:=ZBasis(PI);
    Q,q:=QuotientVS(I,P*I,P);
    queue:={@ Q @};
    output:={@  @};
    done:={@ @};
    while #queue gt 0 do
        pot_new:=&join[MaximalSubmodules(elt) : elt in queue ];
        pot_new_lifts:=[ Ideal(S,[ (Q!b)@@q : b in Basis(W) ] cat zbPI) : W in pot_new]; 
        pot_new:={@ pot_new[i] : i in [1..#pot_new] | O!!pot_new_lifts[i] eq IO @}; //we keep only the ones with trivial extension
        output join:={@ pot_new[i] : i in [1..#pot_new] | not pot_new[i] in done and MultiplicatorRing(pot_new_lift[i]) eq S @};
        done join:=queue;
        // Note: if O!!pot_new_lift[i] is not IO, then all the subvector spaces of pot_new[i] will also not have trivial extension IO.
        // Hence we don't need to continue the recursion on pot_new_lift[i].
        queue := pot_new diff done; 
    end while;
    return output;
end intrinsic;

intrinsic WKICM_bar(S::AlgEtQOrd : Method:="Auto") -> SeqEnum
{Returns all the weak eq classes I, such that (I:I)=S. The VarArg Method (default "Auto") determines if we should use the "IntermediateIdeals" routine or the "LowIndexProcess", which is potentially much slower but more memory efficient.}
    if not assigned S`WKICM_bar then
        if IsGorenstein(S) then
            S`WKICM_bar:=[OneIdeal(S)];
        elif CohenMacaulayType(S) eq 2 then
            S`WKICM_bar:=wkicm_bar_CM_type2(S,NonGorensteinPrimes(S));
        else
            require Method in {"Auto","LowIndexProcess","IntermediateIdeals","IntermediateIdealsVSWithTrivialExtensionAndPrescribedMultiplicatorRing"} : "The VarArg parameter Method is assigned to a not avaialble value";
            vprintf WkClasses,2:"Order of CohenMacaulayType = %o\n",CohenMacaulayType(S);
            // general case
            A:=Algebra(S);
            degA:=Dimension(A);
            seqWk_bar:=[];
            St:=TraceDualIdeal(S);
            if Method in {"Auto","IntermediateIdealsVSWithTrivialExtensionAndPrescribedMultiplicatorRing"} then
                pp:=PrimesAbove(Conductor(S));
                oo:=FindOverOrders(S);
                mult_pp:=[ oo[Index(oo,MultiplicatorRing(P))] : P in pp ];
                assert forall{T : T in mult_pp | assigned T`WKICM_bar};
                num_sub_vect_sp:=function(n,q)
                // q a prime power. Returns the number of F_q-subvector spaces of F_q^n
                // CHECK ME
                    return 2+&+[ (1-q^n)*(1-q) div (1-q^k)*(1-q^(n-k)) : k in [1..n-1]];
                end function;
                sub_vs_T:=[];
                for iP->P in [1..#pp] do
                  wkT:=WKICM_bar(mult_pp[iP]);
                  q:=Inde(S,P);
                  dimsT:=[ Ilog(q,Index(J,(mult_pp[iP]!!P)*J)) : J in wkT ];
                  Append(~sub_vs_T,&+[num_sub_vect_sp(q,d) : d in dimsT]);
                end for;
                _,iP:=Min(sub_vs_T); // this is the T which make us compute less vector spaces.
                P:=pp[iP];
                T:=mult_pp[iP];
                wkT:=WKICM_bar(mult_pp[iP]);
                cands:=[ IntermediateIdealsVSWithTrivialExtensionAndPrescribedMultiplicatorRing(J,P,T), J in wkT ];
                for I in cand do
                  if not exists{J : J in seqWk_bar | IsWeakEquivalent(I,J)} then
                    ZBasisLLL(I);
                    Append(~seqWk_bar,I);
                  end if;
                end for;
            else
              T:=&meet([ T : T in FindOverOrders(S) | IsInvertible(T !! St) ]);
              //this construction of T is conjectural, hence the next assert. If the assert fails, please report it.
              assert IsInvertible(T !! St);

              T1:=S!!OneIdeal(T);
              ff:=ColonIdeal(S,T1); //the relative conductor (S:T)
              Q,q:=Quotient(T1,ff);
              QispGroup:=IspGroup(Q);
              if (Method eq "Auto" and not QispGroup) or (Method eq "IntermediateIdeals") then
                  // // NEW CODE based on MinimalIntermediateIdeals
                  // queue:={@ ff @};
                  // seqWk_bar:=[ OneIdeal(S) , St  ];
                  // done:={@ @};
                  // while #queue gt 0 do
                  //     pot_new:=&join[MinimalIntermediateIdeals(T1,elt) : elt in queue ];
                  //     for I in pot_new do
                  //         if MultiplicatorRing(I) eq S and not exists{ J : J in seqWk_bar | IsWeakEquivalent(I,J) } then
                  //             Append(~seqWk_bar,I);
                  //         end if;
                  //     end for;
                  //     done join:=queue;
                  //     queue := pot_new diff done;
                  // end while;
                  
                  // // NEW CODE based on IntermediateIdeals
                  // idls:=IntermediateIdeals(T1,ff);
                  // for I in idls do
                  //     if MultiplicatorRing(I) eq S and not exists{ J : J in seqWk_bar | IsWeakEquivalent(I,J) } then
                  //         Append(~seqWk_bar,I);
                  //     end if;
                  // end for;
                    
                  // // NEW CODE based on IntermediateIdealsWithPrescribedMultiplicatorRing
                  // vprintf WkClasses,2:"Using code with IntermediateIdealsWithPrescribedMultiplicatorRing: T/(S:T) = %o\n",Q;
                  // idls:=IntermediateIdealsWithPrescribedMultiplicatorRing(T1,ff);
                  // NEW CODE based on IntermediateIdealsWithTrivialExtensionAndPrescribedMultiplicatorRing
                  vprintf WkClasses,2:"Using code with IntermediateIdealsWithTrivialExtensionAndPrescribedMultiplicatorRing: T/(S:T) = %o\n",Q;
                  idls:=IntermediateIdealsWithTrivialExtensionAndPrescribedMultiplicatorRing(T1,ff,T);
                  for I in idls do
                      if not exists{ J : J in seqWk_bar | IsWeakEquivalent(I,J) } then
                          ZBasisLLL(I);
                          Append(~seqWk_bar,I);
                      end if;
                  end for;
              elif (Method eq "Auto" and QispGroup) or (Method eq "LowIndexProcess") then
                  // In this case every subgroup is a module. This seems to be faster
                  //OLD CODE
                  gens_ff_over_S:=Generators(ff);
                  vprintf WkClasses,2:"Using code with LowIndexProcess: T/(S:T) = %o\n",Q;
                  QP,f:=FPGroup(Q); //f:QP->Q
                  subg:=LowIndexProcess(QP,<1,#QP>);
                  while not IsEmpty(subg) do
                      H := ExtractGroup(subg);
                      NextSubgroup(~subg);
                      I:=Ideal(S, [ (f(QP!x))@@q : x in Generators(H) ] cat gens_ff_over_S);
                      if not I in seqWk_bar and 
                          MultiplicatorRing(I) eq S and 
                          not exists{J : J in seqWk_bar | IsWeakEquivalent(I,J)} then 
                              ZBasisLLL(I);
                              Append(~seqWk_bar,I);
                      end if;
                  end while;
              end if;
            end if;
            S`WKICM_bar:=seqWk_bar;
        end if;

    end if;
    return S`WKICM_bar;
end intrinsic;

intrinsic WKICM(E::AlgEtQOrd : Method:="Auto")->SeqEnum
{new}
    require Method in {"Auto","LowIndexProcess","IntermediateIdeals","IntermediateIdealsVSWithTrivialExtensionAndPrescribedMultiplicatorRing"} : "The VarArg parameter Method is assigned to a not avaialble value";
    if not assigned E`WKICM then
        wk:=[];
        oo:=FindOverOrders(E : populateoo_in_oo:=true);
        edges:=[];
        for j in [1..#oo] do
            for i in [Index(oo,T) : T in MinimalOverOrders(oo[j])] do
            // oo[j] subset oo[i] 
                Append(~edges,[i,j]); //the edge is a reversed inclusion
            end for;
        end for;
        D:=Digraph<#oo|edges>;
        vv:=VertexSet(D);
        max_dist:=Distance(vv.Index(oo,MaximalOrder(A)),vv.Index(oo,E));
        for d in [0..max_dist] do
            oo_d:=Sphere(vv.Index(oo,MaximalOrder(A)),d);
            oo_d:={@ oo[Index(v)] : v in oo_d @};
            for T in oo_d do
              wk cat:=[(E!!I) : I in WKICM_bar(T : Method:=Method)];
            end for;
        end for;
        for I in wk do
            ZBasisLLL(I);
        end for;
        E`WKICM:=wk;
    end if;
    return E`WKICM;
end intrinsic;


/* TESTS
    screen -r test_new_wk_icm

    quit;
    cd ~/packages_github/AlgEt/
    git pull; sleep 1;
    magma

    "NEW method";
    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/dev/new_wk_icm.m");

  	_<x>:=PolynomialRing(Integers());
    f:=x^8+16; 
    AttachSpec("~/packages_github/AlgEt/spec");
    A:=EtaleAlgebra(f);
    R:=EquationOrder(A);
    time _:=FindOverOrders(R);
    time #WKICM(R);
    quit;

    magma

    AttachSpec("~/packages_github/AlgEt/spec");
    "OLD method";
  	_<x>:=PolynomialRing(Integers());
    f:=x^8+16; 
    AttachSpec("~/packages_github/AlgEt/spec");
    A:=EtaleAlgebra(f);
    R:=EquationOrder(A);
    time _:=FindOverOrders(R);
    time #WKICM(R);
    quit;

*/
