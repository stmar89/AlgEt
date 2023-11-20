/* vim: set syntax=magma :*/

//freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, stefano.marseglia89@gmail.com
// https://stmar89.github.io/index.html
/////////////////////////////////////////////////////

// This is a variation of rec_wk_icm.m, started while finishing up the paper on July 2023.
// It contains mainly 2 improvements: 
// DONE: in rec_wk_icm_v2.m 1) if R has more than one singular prime ideals, WKICM_bar(R) splits the computation over R+P^nPO, as for WKICM. (at the time of writing rec_wk_icm.m I wasn't sure that this was true.)
// 2) DONE: if R has one singular prime, in the recursion using the extension RcT, we use keep track of the orbits of (T_P)*/(R_P)* rather than lifting the sub-spaces and testing for weak equivalence.
// in the *_v4.m we have included the improvement in the WKICM_bar(S) that if S^tT is invertible, then all fractional ideals in WKICM_bar(S) will extend to the invertible wk. class of T.

//------------
// TODO:
// - remove Method
// - IMPORTANT: once this is done I will have to update also the ICM!!!!
// - change the verbose variable when incorporating in package
declare verbose WKICM, 3;
declare verbose WKICM_bar, 2;
// - when incorporating in package, the next two lines need to be removed (replaced by the ones commented out below)
import "../AlgEtQ/LowCohenMacaulayType.m" : wkicm_bar_CM_type2;
import "../AlgEtQ/PicardGroup.m" : residue_class_ring_unit_subgroup_generators;
//------------

//declare verbose WkClasses, 3;
//declare attributes AlgEtQOrd: WKICM,WKICM_bar;
//import "LowCohenMacaulayType.m" : wkicm_bar_CM_type2;
//import "PicardGroup.m" : residue_class_ring_unit_subgroup_generators;


join_ass_arr:=function(A,B)
// returns A join B, where A[i], B[i] are of type SeqEnum
    for k in Keys(B) do
        if not IsDefined(A,k) then
            A[k]:=B[k];
        else
            for b in B[k] do
                Include(~A[k],b);
            end for;
        end if;
    end for;
    return A;
end function;

diff_ass_arr:=function(A,B,stop_at_dim_2)
// returns A diff B, where A[i], B[i] are of type SeqEnum
// if stop_at_dim_2 then A[2] is removed
    if stop_at_dim_2 then
        Remove(~B,2);
        Remove(~A,2);
    end if;
    for k in Keys(B) do
        if IsDefined(A,k) then
            for b in B[k] do
                Exclude(~A[k],b);
            end for;
        end if;
    end for;
    for k in Keys(A) do
        if #A[k] eq 0 then
            Remove(~A,k);
        end if;
    end for;
    return A;
end function;

units_T_P_mod_S_P:=function(T,S,P)
    // Let S be an order with a unique singular prime ideal P, and T = (P:P)
    // we compute U:=transversal in T of U:=(T_P)*/(S_P)*
    // U = (T/ff)*/(R/ff)*, where ff is the conductor of S in O=MaximalOrder
    // Alternatively, one could use that ff = (S:T) = P, but at the moment there is no implementation to compute A^*, where A=T/P.
    FPS:=Conductor(S);
    O:=MaximalOrder(Algebra(S));
    FPT:=T!!FPS;
    FPO:=O!!FPS;
    uOP,map:=ResidueRingUnits(O,FPO);
    uTP:=sub<uOP|[ g@@map : g in residue_class_ring_unit_subgroup_generators(FPT)]>;
    uSP:=sub<uTP|[uTP!(g@@map) : g in residue_class_ring_unit_subgroup_generators(FPS)]>;
    U:=[ map(uOP!t) : t in Transversal(uTP,uSP) ];
    return U;
end function;

wkicm_bar_with_P_P:=function(I,P)
// Let S be an order, P a prime of S, and I a fractional (P:P)-ideal. 
// The function returns all fractional S-ideals J such that  P*I c J c I, (J:J)=S, and J(P:P)=I, up to weak equivalence.
// We use AssociativeArrays indexed by the dimension at P to make the search faster.
// Also, if I is invertible in (P:P) then we do not consider in the recursion the 1 dimensional vector spaces, 
// since they will all be weakly equivalent to S (which is added manually to the output).
    S:=Order(P);
    T:=Order(I);
    stop_at_dim_2:=false;
    if IsInvertible(I) then
        I:=OneIdeal(T);
        stop_at_dim_2:=true;
    end if;
    TP:=T!!P;
    PI:=TP*I;
    IS:=S!!I;
    zbPI:=ZBasis(PI);
    Q,q:=QuotientVS(IS,P*IS,P); // q:I->I/PI=Q
    if I eq OneIdeal(T) then
        maximal_sub_T_mod:=PrimesAbove(TP); // if I=T then the maximal T-modules P c M c T are precisely
                                            // the primes of T above P. This function is a bit faster than
                                            // the next one, which works also when I ne T.
    else
        maximal_sub_T_mod:=MaximalIntermediateIdeals(I,PI);
    end if;
    maximal_sub_T_mod:=[ sub<Q | [q(z) : z in ZBasis(M)]> : M in maximal_sub_T_mod ] ; // maximal sub-T-modules m of Q,
                                                                                       // whose lift M=q^-1(m) c I
                                                                                       // satisfies PI c M c I
    U:=units_T_P_mod_S_P(T,S,P);

    queue:=AssociativeArray();
    queue[Dimension(Q)]:=[Q];
    output:=AssociativeArray();
    output_vs:=AssociativeArray(); //will contain the orbits
    if MultiplicatorRing(I) eq S then
        output[Dimension(Q)]:=[IS];
    end if;
    if stop_at_dim_2 then
        output[1]:=[OneIdeal(S)];
    end if;
    done:=AssociativeArray();
    while #queue gt 0 do //empty keys are removed
        // We want only the sub-S/P-vector spaces W of Q such that W.T=Q (where ".T" denotes the natural action).
        // This is equivalent to having JT=I, where J=q^-1(W) c I.
        // We use that: JT=I if and only if W is not contained in any M for M in maximal_sub_T_mod.
        // Also, if W does not have trivial extension, then all its subspace will not have trivial extension.
        // Hence we do not add such W to the queue.
        pot_new:=AssociativeArray();
        for W in &cat[ &cat[ MaximalSubmodules(elt) : elt in queue[k_queue] ] : k_queue in Keys(queue) ] do
            dimW:=Dimension(W);
            if not IsDefined(pot_new,dimW) then
                pot_new[dimW]:=[];
            end if;
            if not W in pot_new[dimW] and not exists{ M : M in maximal_sub_T_mod | W subset M } then
                Append(~pot_new[dimW],W);
                J:=Ideal(S,[ (Q!b)@@q : b in Basis(W) ] cat zbPI);
                if MultiplicatorRing(J) eq S then
                    if not IsDefined(output_vs,dimW) or not W in output_vs[dimW] then
                        // something new! we compute the orbit
                        zbJ:=Generators(J);
                        orb_vs:=[ sub<Q | [q(u*g) : g in zbJ]> : u in U ];
                        // we add the orbit to the output_vs, and J to the output
                        for i in [1..#orb_vs] do
                            dimJJ:=Dimension(orb_vs[i]);  
                            if not IsDefined(output_vs,dimJJ) then
                                output_vs[dimJJ]:=[orb_vs[i]];
                            else
                                Append(~output_vs[dimJJ],orb_vs[i]);
                            end if;
                        end for;
                        if not IsDefined(output,dimW) then
                            output[dimW]:=[J];
                        else
                            Append(~output[dimW],J);
                        end if;
                    end if;
                end if;
            end if;
        end for;
        done:=join_ass_arr(done,queue); // done join:=queue;
        queue:=diff_ass_arr(pot_new,done,stop_at_dim_2); // queue := pot_new diff done;
    end while;
    output:=&cat[ output[k] : k in Keys(output) ];
    return output;
end function;

intrinsic WKICM_bar(S::AlgEtQOrd : Method:="Auto") -> SeqEnum[AlgEtQIdl]
{Returns all the weak eq classes I, such that (I:I)=S. The VarArg Method (default "Auto") determines if we should use the "IntermediateIdeals" routine or the "LowIndexProcess", which is potentially much slower but more memory efficient.}
    if not assigned S`WKICM_bar then
        if IsGorenstein(S) then
            vprintf WKICM_bar,2 : " Gorenstein case\n";
            S`WKICM_bar:=[OneIdeal(S)];
        elif CohenMacaulayType(S) eq 2 then
            vprintf WKICM_bar,2 : " Cohen Macaulay type 2 case\n";
            S`WKICM_bar:=wkicm_bar_CM_type2(S,NonGorensteinPrimes(S));
        else
            vprintf WkClasses,2:"Order of CohenMacaulayType = %o\n",CohenMacaulayType(S);
            // general case
            seqWk_bar:=[];
            vprint WKICM_bar,2 : "Using new method";
            pp:=SingularPrimes(S);
            if #pp eq 1 then
                vprint WKICM_bar,2 : "WKICM_bar: type >2, only 1 singular prime";
                P:=pp[1];
                T:=MultiplicatorRing(P);
                // if S^tT is invertible, then I don't need all WKICM_bar(T), but only the invertible class.
                if IsInvertible(T!!TraceDualIdeal(S)) then
                    wkT:=[ OneIdeal(T) ];
                else
                    wkT:=$$(T);
                end if;
                for J in wkT do
                    seqWk_bar_J:=wkicm_bar_with_P_P(J,P);
                    for I in seqWk_bar_J do
                        ZBasisLLL(I);
                    end for;
                    Append(~seqWk_bar,seqWk_bar_J);
                end for;
                vprintf WKICM_bar,2 : "sizes of seqWk_bar_J = %o\n",[#x : x in seqWk_bar];
                seqWk_bar:=&cat(seqWk_bar);
                S`WKICM_bar:=seqWk_bar;
            else
                vprint WKICM_bar,2 : "WKICM_bar: type >2, more than 1 singular prime";
                O:=MaximalOrder(Algebra(S));
                wk_pp:=[];
                // We have W_bar(S) = \prod_P W_P_bar(S), where:
                //     -  P runs over the singular primes of S, 
                //     - W_bar(S) is the monoid of weak equivalence classes with mult ring S, and
                //     - W_P_bar(S) is monoid of local P-equivalence classes with mult ring equal to S at P.
                // We also have W_P_bar(S) = W_bar(S+P^kO), for any k big enough (eg. k ge v_p([O:S]),
                //       where p is the rational prime of P).
                ps:=[]; // rational primes of the pps
                for iP->P in pp do
                    t0:=Cputime();
                    vprintf WKICM_bar,2 : "We start the local computation at the %o-th singular prime.",iP;
                    wk_P:=[];
                    _,p:=IsPrimePower(Index(S,P));
                    Append(~ps,p);
                    k:=Valuation(Index(O,S),p);
                    SP:=Order(ZBasis(S) cat ZBasis(O!!P^k));
                    wk_P_bar:=$$(SP); // recursion
                    // wk_P = W_P(S) = W(S+P^k0O) in the notation above
                    wk_P_bar:=[ S !! I : I in wk_P_bar];
                    vprintf WKICM_bar,2 : "...Done in %o secs.\n",Cputime(t0);
                    Append(~wk_pp,wk_P_bar);
                end for;

                ///////////////////////////////
                // We now reconstruct W_bar(S) using the P-local parts W_P_bar(S) = W_bar(S+P^kO)
                //////////////////////////////
                wk_pp_idls:=[];
                pp_pows:=[];
                t1:=Cputime();
                vprintf WKICM_bar,2 : "We make all the local parts integral\n";
                for ip->wk in wk_pp do
                    wk_exps:=[];
                    wk_idls:=[];
                    for i in [1..#wk] do
                        I:=wk[i];
                        if not IsIntegral(I) then
                            I:=SmallRepresentative(I); // I c E with small norm
                        end if;
                        k:=Valuation(Index(S,I),ps[ip]);
                        Append(~wk_exps,k);
                        Append(~wk_idls,I);
                    end for;
                    k_ip:=Max(wk_exps);
                    Pk_ip:=pp[ip]^k_ip; // for every local representative I at pp[ip] we have that Pk_ip c I (locally)
                    ZBasisLLL(Pk_ip);
                    Append(~pp_pows,Pk_ip);
                    Append(~wk_pp_idls,wk_idls);
                end for;
                vprintf WKICM_bar,2 : "...Done in %o secs.\n",Cputime(t1);
                    
                n:=#pp;
                t0:=Cputime();
                vprintf WKICM_bar,2 : "We compute the \prod_{j \\ne i} P_j^k_j\n";
                prod_j_ne_i:=[ ];
                for i in [1..n] do
                    prod:=&*[ pp_pows[j] : j in [1..n] | j ne i ];
                    ZBasisLLL(prod);
                    Append(~prod_j_ne_i,prod);
                end for;
                vprintf WKICM_bar,2 : "\t...Done in %o secs.\n",Cputime(t0);

                t0:=Cputime();
                vprintf WKICM_bar,2 : "We modify each entry of the cartesian product\n";
                for ip in [1..n] do
                    for i in [1..#wk_pp_idls[ip]] do
                        I:=(wk_pp_idls[ip][i]+pp_pows[ip])*prod_j_ne_i[ip];
                        ZBasisLLL(I);
                        wk_pp_idls[ip][i]:=I;
                    end for;
                end for;
                vprintf WKICM_bar,2 : "\t...Done in %o secs.\n",Cputime(t0);

                t0:=Cputime();
                tot:=&*[#x : x in wk_pp_idls]; perc_old:=0; iI:=0;
                wk_pp_idls:=CartesianProduct(wk_pp_idls);
                vprintf WKICM_bar,2 : "We start patching together the local parts\n";
                wk:=[];
                for I_Ps in wk_pp_idls do
                    if GetVerbose("WKICM_bar") ge 3 then
                        iI +:=100; perc:=Truncate(iI/tot); 
                        if perc gt perc_old then perc_old:=perc; printf "\t%o%% in %o secs\n",perc,Cputime(t0); end if;
                    end if;
                    J:=&+[ I_Ps[ip] : ip in [1..n] ];
                    // J satisfies: J = I_Ps[ip] locally at pp[ip] for every ip.
                    assert2 forall{ ip : ip in [1..n] | 
                                                    (J+I_Ps[ip]) eq I_Ps[ip]+pp[ip]*(J+I_Ps[ip]) and 
                                                    (J+I_Ps[ip]) eq J+pp[ip]*(J+I_Ps[ip])};
                    Append(~wk,J);
                end for;
                vprintf WKICM_bar,2 : "\t...Done in %o secs.\n",Cputime(t0);

                t0:=Cputime();
                vprintf WKICM_bar,2 : "We LLL all the ZBasis\n";
                for I in wk do
                    ZBasisLLL(I);
                end for;
                vprintf WKICM_bar,2 : "\t...Done in %o secs\n",Cputime(t0);

                t0:=Cputime();
                vprintf WKICM_bar,2 : "Checking assert2 ...\n";
                // asserts
                assert2 forall{ J : J in wk | Order(J) eq S };
                assert2 forall{ J : J in wk | not exists{I : I in wk | I ne J and IsWeakEquivalent(I,J) }  };
                vprintf WKICM_bar,2 : "\t...Done in %o secs\n",Cputime(t0);

                S`WKICM_bar:=wk;
            end if;
        end if;
    end if;
    // populate MultiplicatorRing, if not done before
    for i in [1..#S`WKICM_bar] do
        I:=S`WKICM_bar[i];
        if not assigned I`MultiplicatorRing then
            I`MultiplicatorRing:=S;
        end if;
    end for;
    return S`WKICM_bar;
end intrinsic;

intrinsic WeakEquivalenceClassesWithPrescribedMultiplicatorRing(S::AlgEtQOrd : Method:="Auto") -> SeqEnum[AlgEtQIdl]
{}
    return WKICM_bar(S : Method:=Method);
end intrinsic;

intrinsic WKICM(E::AlgEtQOrd : Method:="Auto")->SeqEnum[AlgEtQIdl]
{}
    require Method in {"Auto","LowIndexProcess","IntermediateIdeals","IntermediateIdealsVSWithTrivialExtensionAndPrescribedMultiplicatorRing"} : "The VarArg parameter Method is assigned to a not available value";
    if not assigned E`WKICM then
        pp:=SingularPrimes(E);
        if #pp eq 0 then
            E`WKICM:=[OneIdeal(E)]; // the order is maximal
            E`WKICM_bar:=[OneIdeal(E)]; // the order is maximal
        elif #pp eq 1 then // 1 singular prime: go for the recursion
            vprintf WKICM,2 : "There is only one singular prime.\n";
            // W(E) = union_{ S in MinimalOverOrders(E)} W(S) disjoint_union WKICM_bar(E);
            // The first union should be made in a smart way!
            min:=MinimalOverOrders(E);
            for S in min do
                _:=$$(S); // recursion to populate
            end for;
            // We want to put in oo_E the strict overorders of E, from the overorders of the minimal overorders.
            // There will be a lot of overlap. So we try to do it in a smart way.
            t0:=Cputime();
            vprintf WKICM,2 : "\tWe merge the overorders of the minimal overorders.";
            oo_E:=AssociativeArray();
            for S in min do
                for T in OverOrders(S) do
                    indT:=Index(T);
                    if not IsDefined(oo_E,indT) then
                        oo_E[indT]:={@ T @};
                    else
                        Include(~oo_E[indT],T);
                    end if;
                end for;
            end for;
            oo_E:=&cat[ Setseq(oo_E[k]) : k in Keys(oo_E) ];
            vprintf WKICM,2 : "...Done in %o secs.\n",Cputime(t0);
            assert2 forall{ S : S in oo_E | assigned S`WKICM_bar };
            wk_mins:=&cat[ [ E !! I : I in WKICM_bar(S)] : S in oo_E ];
            E`OverOrders:=[ E ] cat oo_E;
            E`WKICM:=wk_mins cat WKICM_bar(E); // here I need to use WKICM_bar(P:P) which is already computed since
                                               // each for each S in min we have that S c (P:P)
        else //more then 1 singular prime: divide and patch
            vprintf WKICM,2 : "There are %o singular primes.\n",#pp;
            O:=MaximalOrder(Algebra(E));
            wk_pp:=[];
            // We have W(E) = \prod_P W_P(E), where:
            //     -  P runs over the singular primes of E, 
            //     - W(E) is the monoid of weak equivalence classes, and 
            //     - W_P(E) is monoid of local P-equivalence classes.
            // We also have W_P(E) = W(E+P^kO), for any k big enough (eg. k ge v_p([O:E]),
            //       where p is the rational prime of P).
            ps:=[]; // rational primes of the pps
            for iP->P in pp do
                t0:=Cputime();
                vprintf WKICM,2 : "We start the local computation at the %o-th singular prime.",iP;
                wk_P:=[];
                _,p:=IsPrimePower(Index(E,P));
                Append(~ps,p);
                k:=Valuation(Index(O,E),p);
                EP:=Order(ZBasis(E) cat ZBasis(O!!P^k));
                wk_P:=$$(EP); // recursion
                // wk_P = W_P(E) = W(E+P^k0O) in the notation above
                wk_P:=[ E !! I : I in wk_P];
                vprintf WKICM,2 : "...Done in %o secs.\n",Cputime(t0);
                Append(~wk_pp,wk_P);
            end for;

            ///////////////////////////////
            // We now reconstruct W(E) using the P-local parts W_P(E) = W(E+P^kO)
            //////////////////////////////
            wk_pp_idls:=[];
            pp_pows:=[];
            t1:=Cputime();
            vprintf WKICM,2 : "We make all the local parts integral\n";
            for ip->wk in wk_pp do
                wk_exps:=[];
                wk_idls:=[];
                for i in [1..#wk] do
                    I:=wk[i];
                    if not IsIntegral(I) then
                        I:=SmallRepresentative(I); // I c E with small norm
                    end if;
                    k:=Valuation(Index(E,I),ps[ip]);
                    Append(~wk_exps,k);
                    Append(~wk_idls,I);
                end for;
                k_ip:=Max(wk_exps);
                Pk_ip:=pp[ip]^k_ip; // for every local representative I at pp[ip] we have that Pk_ip c I (locally)
                ZBasisLLL(Pk_ip);
                Append(~pp_pows,Pk_ip);
                Append(~wk_pp_idls,wk_idls);
            end for;
            vprintf WKICM,2 : "...Done in %o secs.\n",Cputime(t1);
                
            n:=#pp;
            t0:=Cputime();
            vprintf WKICM,2 : "We compute the \prod_{j \\ne i} P_j^k_j\n";
            prod_j_ne_i:=[ ];
            for i in [1..n] do
                prod:=&*[ pp_pows[j] : j in [1..n] | j ne i ];
                ZBasisLLL(prod);
                Append(~prod_j_ne_i,prod);
            end for;
            vprintf WKICM,2 : "\t...Done in %o secs.\n",Cputime(t0);

            t0:=Cputime();
            vprintf WKICM,2 : "We modify each entry of the cartesian product\n";
            for ip in [1..n] do
                for i in [1..#wk_pp_idls[ip]] do
                    I:=(wk_pp_idls[ip][i]+pp_pows[ip])*prod_j_ne_i[ip];
                    ZBasisLLL(I);
                    wk_pp_idls[ip][i]:=I;
                end for;
            end for;
            vprintf WKICM,2 : "\t...Done in %o secs.\n",Cputime(t0);

            t0:=Cputime();
            tot:=&*[#x : x in wk_pp_idls]; perc_old:=0; iI:=0;
            wk_pp_idls:=CartesianProduct(wk_pp_idls);
            vprintf WKICM,2 : "We start patching together the local parts\n";
            wk:=[];
            for I_Ps in wk_pp_idls do
                if GetVerbose("WKICM") ge 3 then
                    iI +:=100; perc:=Truncate(iI/tot); 
                    if perc gt perc_old then perc_old:=perc; printf "\t%o%% in %o secs\n",perc,Cputime(t0); end if;
                end if;
                J:=&+[ I_Ps[ip] : ip in [1..n] ];
                // J satisfies: J = I_Ps[ip] locally at pp[ip] for every ip.
                assert2 forall{ ip : ip in [1..n] | 
                                                (J+I_Ps[ip]) eq I_Ps[ip]+pp[ip]*(J+I_Ps[ip]) and 
                                                (J+I_Ps[ip]) eq J+pp[ip]*(J+I_Ps[ip])};
                Append(~wk,J);
            end for;
            vprintf WKICM,2 : "\t...Done in %o secs.\n",Cputime(t0);

            t0:=Cputime();
            vprintf WKICM,2 : "We LLL all the ZBasis\n";
            for I in wk do
                ZBasisLLL(I);
            end for;
            vprintf WKICM,2 : "\t...Done in %o secs\n",Cputime(t0);

            t0:=Cputime();
            vprintf WKICM,2 : "Checking assert2 ...\n";
            // asserts
            assert2 forall{ J : J in wk | Order(J) eq E };
            assert2 forall{ J : J in wk | not exists{I : I in wk | I ne J and IsWeakEquivalent(I,J) }  };
            vprintf WKICM,2 : "\t...Done in %o secs\n",Cputime(t0);

            E`WKICM:=wk;

            // We populate the attributes E`OverOrders and S`WKICM_bar for each overorder.
            // This is needed for the recursion step.
            // Note that some of the overorders have already the vararg WKICM_bar populated, because, for 
            // example, there were coming from the recursive-1-single-singular-prime part of the code.
            // In this case, we don't want to add it a second time. Hence we need to keep track of this info.
            // We do it by using nested AossicativeArrays
            t0:=Cputime();
            vprintf WKICM,2 : "We start populating the varargs OverOrders and WKICM_bar.\n";
            if not assigned E`OverOrders then
                oo:=AssociativeArray();
                for I in wk do
                    S:=MultiplicatorRing(I);
                    ind:=Index(S);
                    if not IsDefined(oo,ind) then
                        oo[ind]:=AssociativeArray();
                        oo[ind][true]:={@ @};
                        oo[ind][false]:=AssociativeArray();
                        oo[ind][false]["orders"]:=[];
                        oo[ind][false]["wkicm_bars"]:=[];
                    end if;
                    if assigned S`WKICM_bar then
                        Include(~oo[ind][true],S);
                    else
                        indS:=Index(oo[ind][false]["orders"],S);
                        if indS eq 0 then //this is the first time we see S.
                            Append(~oo[ind][false]["orders"],S);
                            Append(~oo[ind][false]["wkicm_bars"],[S!!I]);
                        else
                            Append(~oo[ind][false]["wkicm_bars"][indS],S!!I);
                        end if;
                    end if;
                end for;
                oo_E:=[];
                for ind in Keys(oo) do
                    oo_E cat:=Setseq(oo[ind][true]);
                    for iS in [1..#oo[ind][false]["orders"]] do
                        S:=oo[ind][false]["orders"][iS];
                        S`WKICM_bar:=oo[ind][false]["wkicm_bars"][iS];
                        Append(~oo_E,S);
                    end for;
                end for;
                E`OverOrders:=oo_E;
            end if;
            vprintf WKICM,2 : "\t...Done in %o secs.\n",Cputime(t0);
            assert2 &+[ #S`WKICM_bar : S in OverOrders(E) ] eq #wk;
        end if;
    end if;

    return E`WKICM;
end intrinsic;

/* TESTS
    see rec_wk_icm_v2_timings.m
*/

