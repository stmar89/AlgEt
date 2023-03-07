/* vim: set syntax=magma :*/

//freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////


//------------
// TODO:
// - IMPORTANT: once this is done I will have to update also the ICM!!!!
// - change the verbose variable when incorporating in package
declare verbose WKICM, 3;
declare verbose WKICM_bar, 2;
// - when incorporating in package, the next line needs to be removed (replaced by the one commented out below)
import "../AlgEtQ/LowCohenMacaulayType.m" : wkicm_bar_CM_type2;
//------------

//declare verbose WkClasses, 3;
//declare attributes AlgEtQOrd: WKICM,WKICM_bar;
//import "LowCohenMacaulayType.m" : wkicm_bar_CM_type2;

/* not used anymore
intrinsic IntermediateIdealsVSWithTrivialExtensionAndPrescribedMultiplicatorRing(I::AlgEtQIdl,P::AlgEtQIdl)->SetIndx[AlgEtQIdl]
{Let S be an order, P a prime of S, and I a fractional (P:P)-ideal. The intrinsic returns all fractional S-ideals J such that  P*I c J c I, (J:J)=S, and J(P:P)=I.}
    require IsPrime(P) : "The ideal P must be prime";
    S:=Order(P);
    T:=Order(I);
    require T eq MultiplicatorRing(P): "The ideal I must be over (P:P)";
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
    queue:={@ Q @};
    output:={@  @};
    if MultiplicatorRing(I) eq S then
        Include(~output, IS);
    end if;
    done:={@ @};
    while #queue gt 0 do
        pot_new:={@ W : W in &cat[ MaximalSubmodules(elt) : elt in queue ]@};
        pot_new:={@ W : W in pot_new | not exists{ M : M in maximal_sub_T_mod | W subset M} @};
            // We want only the sub-S/P-vector spaces W of Q such that W.T=Q (where ".T" denotes the natural action).
            // This is equivalent to having JT=I, where J=q^-1(W) c I.
            // We use that: JT=I if and only if W is not contained in any M for M in maximal_sub_T_mod.
            // Also, if W does not have trivial extension, then all its subspace will not have trivial extension.
            // Hence we do not add such W to the queue.
        pot_new_lifts:=[ Ideal(S,[ (Q!b)@@q : b in Basis(W) ] cat zbPI) : W in pot_new]; 
        output join:={@ pot_new_lifts[i] : i in [1..#pot_new] | MultiplicatorRing(pot_new_lifts[i]) eq S @};
            // one the lifts with multiplicator ring S
        done join:=queue;
        queue := pot_new diff done; 
    end while;
    return output;
end intrinsic;
*/

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
    queue:=AssociativeArray();
    queue[Dimension(Q)]:=[Q];
    output:=AssociativeArray();
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
                I:=Ideal(S,[ (Q!b)@@q : b in Basis(W) ] cat zbPI);
                if MultiplicatorRing(I) eq S then
                    if not IsDefined(output,dimW) then
                        output[dimW]:=[I];
                    elif not exists{ J : J in output[dimW] | IsWeakEquivalent(I,J) } then
                        Append(~output[dimW],I);
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
            require Method in {"Auto","LowIndexProcess","IntermediateIdeals","IntermediateIdealsVSWithTrivialExtensionAndPrescribedMultiplicatorRing"} : "The VarArg parameter Method is assigned to a not avaialble value";
            vprintf WkClasses,2:"Order of CohenMacaulayType = %o\n",CohenMacaulayType(S);
            // general case
            seqWk_bar:=[];
            vprint WKICM_bar,2 : "Using new method";
            pp:=PrimesAbove(Conductor(S));
            mult_pp:=[ MultiplicatorRing(P) : P in pp ];
            assert forall{T : T in mult_pp | assigned T`WKICM_bar};

            if #pp gt 1 then
                num_sub_vect_sp:=function(n,q)
                // q a prime power. Returns the number of F_q-subvector spaces of F_q^n
                    return &+[ GaussianBinomial(n,k,q) : k in [0..n]];
                end function;
                sub_vs_T:=[];
                for iP->P in pp do
                  wkT:=WKICM_bar(mult_pp[iP]);
                  q:=Index(S,P);
                  dimsT:=[ Ilog(q,Index(J,(mult_pp[iP]!!P)*J)) : J in wkT ];
                  Append(~sub_vs_T,&+[num_sub_vect_sp(d,q) : d in dimsT]);
                end for;
                _,iP:=Min(sub_vs_T); // this is the T which make us compute less vector spaces.
            else
                iP:=1;
            end if;
            P:=pp[iP];
            T:=mult_pp[iP];
            wkT:=$$(T);
            for J in wkT do
                /* OLD rec without AssociativeArray
                seqWk_bar_J:=[];
                cands_J:=IntermediateIdealsVSWithTrivialExtensionAndPrescribedMultiplicatorRing(J,P);
                for I in cands_J do
                    if not exists{K : K in seqWk_bar_J | IsWeakEquivalent(I,K)} then
                        ZBasisLLL(I);
                        Append(~seqWk_bar_J,I);
                    end if;
                end for;
                */
                seqWk_bar_J:=wkicm_bar_with_P_P(J,P);
                for I in seqWk_bar_J do
                    ZBasisLLL(I);
                end for;
                Append(~seqWk_bar,seqWk_bar_J);
            end for;
            vprintf WKICM_bar,2 : "sizes of seqWk_bar_J = %o\n",[#x : x in seqWk_bar];
            seqWk_bar:=&cat(seqWk_bar);
            S`WKICM_bar:=seqWk_bar;
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
 
    //SetAssertions(2);

    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/dev/rec_wk_icm.m");
    //SetVerbose("WKICM",2);
    P<x>:=PolynomialRing(Integers());
    f:=x^8+16;
    A:=EtaleAlgebra(f);
    R:=EquationOrder(A);
    // only one singular prime
    t0:=Cputime();
        assert #WKICM(R) eq 173;
    t_curr:=Cputime(t0);
    t_prev_best:=32.7; //22.9
    "Current running time: ",t_curr;
    if Abs((t_curr - t_prev_best)/t_prev_best) gt 0.1 then 
        print "The current timing is different from the previous best known one. UPDATE!"; 
    end if;
    "the size of the output, 173 classes, has been computed using the OLD method, in 189000 seconds";

    // these tests have more than one singular prime. each takes <1 min
    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/dev/rec_wk_icm.m");
    //SetVerbose("WKICM",2);
    P<x>:=PolynomialRing(Integers());
    f:=x^10 - x^9 + 4*x^8 - 6*x^7 + 8*x^6 - 16*x^5 + 16*x^4 - 24*x^3 + 32*x^2 - 16*x + 32;
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    q:=Integers() ! Round(ConstantCoefficient(f)^(2/Degree(f)));
    R:=Order([F,q/F]);
    t0:=Cputime();
        wk:=WKICM(R);
    t_curr:=Cputime(t0);
    "Current running time: ",t_curr;
    t_prev_best:=12.7; //10.6
    if Abs((t_curr - t_prev_best)/t_prev_best) gt 0.1 then 
        print "The current timing is different from the previous best known one. UPDATE!"; 
    end if;

    // Here: R has 3 singular primes, whose local parts require ~60 secs, the patching is <5 secs.
    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/dev/rec_wk_icm.m");
    //SetVerbose("WKICM",2);
    P<x>:=PolynomialRing(Integers());
    f:=x^6 - 3*x^5 + 30*x^4 - 175*x^3 + 750*x^2 - 1875*x + 15625;
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    q:=Integers() ! Round(ConstantCoefficient(f)^(2/Degree(f)));
    R:=Order([F,q/F]);
    t0:=Cputime();
        #WKICM(R);
    t_curr:=Cputime(t0);
    "Current running time: ",t_curr;
    t_prev_best:=57; //31
    if Abs((t_curr - t_prev_best)/t_prev_best) gt 0.1 then 
        print "The current timing is different from the previous best known one. UPDATE!"; 
    end if;

    // Here: R has 5 singular primes, whose local parts require <60 secs, while the patching requires <5 secs.
    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/dev/rec_wk_icm.m");
    //SetVerbose("WKICM",2);
    P<x>:=PolynomialRing(Integers());
    f:=x^6 + 8*x^5 + 50*x^4 + 200*x^3 + 1250*x^2 + 5000*x + 15625;
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    R:=Order([F,25/F]);
    t0:=Cputime();
        wk:=WKICM(R);
    t_curr:=Cputime(t0);
    "Current running time: ",t_curr;
    t_prev_best:=53; //30.9
    if Abs((t_curr - t_prev_best)/t_prev_best) gt 0.1 then 
        print "The current timing is different from the previous best known one. UPDATE!"; 
    end if;

    // VERY big it should not be in this test suit
    // I have improved a lot KnownOrders. I expect that populating OverORders should be much faster now
    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/dev/rec_wk_icm.m");
    //SetVerbose("WKICM",2);
    P<x>:=PolynomialRing(Integers());
    f:=x^8 - 2*x^7 + 7*x^6 - 14*x^5 + 40*x^4 - 56*x^3 + 112*x^2 - 128*x + 256;
    A:=EtaleAlgebra(f);
    R:=EquationOrder(A);
    time assert #WKICM(R) eq 114492; // 114492 in ~519000 second ~ 144h with the OLD method
                                     // Now we have 1.63 h
*/

/*

// old tests

    // this seems to be quite big. I don't know how long it takes
    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/dev/new_wk_icm.m");
    //SetVerbose("new_wk_icm_bar",2);
    P<x>:=PolynomialRing(Integers());
    f:=x^6 + 8*x^5 + 50*x^4 + 200*x^3 + 1250*x^2 + 5000*x + 15625;
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    R:=Order([F,25/F]);
    assert #SingularPrimes(R) eq 5;
    time #WKICM(R);
    quit;

    magma;
    // OLD method
    AttachSpec("~/packages_github/AlgEt/spec");
    //SetVerbose("new_wk_icm_bar",2);
    P<x>:=PolynomialRing(Integers());
    f:=x^6 + 8*x^5 + 50*x^4 + 200*x^3 + 1250*x^2 + 5000*x + 15625;
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    R:=Order([F,25/F]);
    time #WKICM(R); //560 classes 160 secs


    // VERY big
    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/dev/new_wk_icm.m");
    SetVerbose("WKICM",2);
    P<x>:=PolynomialRing(Integers());
    f:=x^8 - 2*x^7 + 7*x^6 - 14*x^5 + 40*x^4 - 56*x^3 + 112*x^2 - 128*x + 256;
    A:=EtaleAlgebra(f);
    R:=EquationOrder(A);
    time assert #WKICM(R) eq 114492; // 114492 in ~519000 second ~ 144h with the old method
    // Computing the local parts takes <8000 secs. But putting them together might take a long time.
    
    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/dev/new_wk_icm.m");
    // to load the huge example:
    time R:=LoadWKICM(Read("~/packages_github/AlgEt/dev/114492_wk_classes_example.txt")); // 3400 secs to load :-)
    #WKICM(R);


    // looping over slow input
    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/dev/new_wk_icm.m");
    //SetVerbose("WKICM",2);
    P<x>:=PolynomialRing(Integers());
    str:=Split(Read("~/packages_github/AlgEt/dev/input_slow_sorted.txt"));
    poly:=[ P!eval(l) : l in str ];
    for i->f in poly do
        A:=EtaleAlgebra(f);
        F:=PrimitiveElement(A);
        q:=Integers() ! Round(ConstantCoefficient(f)^(2/Degree(f)));
        R:=Order([F,q/F]);
        if #SingularPrimes(R) gt 1 then
            printf "f=%o\n",f;
            t0:=Cputime();
            wk:=#WKICM(R);
            printf "\t#wk=%o seconds=%o\n\n",wk,Cputime(t0);
        end if;
    end for;

    // freakking big
    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/dev/new_wk_icm.m");
    //SetVerbose("new_wk_icm_bar",2);
    P<x>:=PolynomialRing(Integers());
    f:=x^8 - 2*x^7 + 7*x^6 - 14*x^5 + 40*x^4 - 56*x^3 + 112*x^2 - 128*x + 256;
    A:=EtaleAlgebra(f);
    R:=EquationOrder(A);
    "Computing overorders ...";
    time oo:=FindOverOrders(R); //14824 seconds, 16320 overorders
    #oo, "overorders found.";
    time wk:=WKICM(R); // 114492i in ~519000 second ~ 144h
    
    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/dev/new_wk_icm.m");
    // to load the huge example:
    time R:=LoadWKICM(Read("~/packages_github/AlgEt/dev/114492_wk_classes_example.txt")); // 3400 secs to load :-)
    #WKICM(R);

    // approach by splitting the computation over the 2 rational singular primes
    wks:=[];
    O:=MaximalOrder(A);
    fac:=Factorization(Index(O,R));
    for ff in fac do
        // Rp:=[ S : S in oo | Index(O,S) eq ff[1]^ff[2] ];
        // assert #Rp eq 1;
        // Rp:=Rp[1];
        Rp:=Order(ZBasis(R) cat ZBasis(ff[1]^ff[2]*O));
        time wkp:=WKICM(Rp); //11000 seconds for p=2. the other one is superfast
        Append(~wks,#wkp);    
    end for;
    assert &*wks eq 114492;

    // approach by splitting the computation over the 3 singular primes.
    // here we consider 3 overorders by looking at primes of R
    // there are two primes P above the most singular rational prime p=2.
    // let's see if it is faster. Yes! the first and third primes have trivial running time
    // the second is around 7000 secs
    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/dev/new_wk_icm.m");
    P<x>:=PolynomialRing(Integers());
    f:=x^8 - 2*x^7 + 7*x^6 - 14*x^5 + 40*x^4 - 56*x^3 + 112*x^2 - 128*x + 256;
    A:=EtaleAlgebra(f);
    O:=MaximalOrder(A);
    R:=EquationOrder(A);
    pp:=SingularPrimes(R);
    wks:=[];
    for P in pp do
        _,p,_:=IsPrimePower(Index(R,P));
        k:=Valuation(Index(O,R),p);
        Rp:=Order(ZBasis(R) cat ZBasis(p^k*O));
        RP:=Order(ZBasis(R) cat ZBasis(O!!P^k));
        Rp subset RP,Index(RP,Rp),Index(O,RP);
        time Append(~wks,#WKICM(RP));
    end for;
    assert &*wks eq 114492;

    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/dev/new_wk_icm.m");
    //SetVerbose("new_wk_icm_bar",2);
    P<x>:=PolynomialRing(Integers());
    f:=x^4+30^3*x^3+30^3*x^2+30^3*x+30^3;
    A:=EtaleAlgebra(f);
    R:=EquationOrder(A);
    wks:=[];
    O:=MaximalOrder(A);
    fac:=Factorization(Index(O,R));
    for ff in fac do
        // Rp:=[ S : S in oo | Index(O,S) eq ff[1]^ff[2] ];
        // assert #Rp eq 1;
        // Rp:=Rp[1];
        Rp:=Order(ZBasis(R) cat ZBasis(ff[1]^ff[2]*O));
        time wkp:=WKICM(Rp);
        Append(~wks,#wkp);    
    end for;
    time #WKICM(R);
    assert &*wks eq #WKICM(R);

    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/dev/new_wk_icm.m");
    //SetVerbose("new_wk_icm_bar",2);
    P<x>:=PolynomialRing(Integers());
    f:=x^6 + 8*x^5 + 50*x^4 + 200*x^3 + 1250*x^2 + 5000*x + 15625;
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    O:=MaximalOrder(A);
    R:=EquationOrder(A); //3200 overorders
    #SingularPrimes(R);
    fac:=Factorization(Index(O,R));
    #fac;
    wks:=[];
    for ff in fac do
        // Rp:=[ S : S in oo | Index(O,S) eq ff[1]^ff[2] ];
        // assert #Rp eq 1;
        // Rp:=Rp[1];
        Rp:=Order(ZBasis(R) cat ZBasis(ff[1]^ff[2]*O));
        time wkp:=WKICM(Rp);
        Append(~wks,#wkp);    
    end for;
    time #WKICM(R);
    assert &*wks eq #WKICM(R);

*/
/*
    // before implementing the prime-per-prime:
    // testing WKICM_LESS_OLD (trivial extenesion with sub-vector-spaces) vs WKICM_OLD (trivial extension with lifts)
    // on diophntus LESS_OLD is faster 8000sec vs 9000 
    // the same seems to hold on gemini 7474 vs 7980
    quit;
    magma;
    // WKICM_LESS_OLD
    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/dev/new_wk_icm.m");
    P<x>:=PolynomialRing(Integers());
    f:=x^8 - 2*x^7 + 7*x^6 - 14*x^5 + 40*x^4 - 56*x^3 + 112*x^2 - 128*x + 256;
    A:=EtaleAlgebra(f);
    O:=MaximalOrder(A);
    R:=EquationOrder(A);
    pp:=SingularPrimes(R);
    wks:=[];
    for P in pp do
        _,p,_:=IsPrimePower(Index(R,P));
        k:=Valuation(Index(O,R),p);
        RP:=Order(ZBasis(R) cat ZBasis(O!!P^k));
        time Append(~wks,#WKICM_LESS_OLD(RP));
    end for;
    assert &*wks eq 114492;

    quit;
    magma;

    // WKICM_OLD
    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/dev/new_wk_icm.m");
    P<x>:=PolynomialRing(Integers());
    f:=x^8 - 2*x^7 + 7*x^6 - 14*x^5 + 40*x^4 - 56*x^3 + 112*x^2 - 128*x + 256;
    A:=EtaleAlgebra(f);
    O:=MaximalOrder(A);
    R:=EquationOrder(A);
    pp:=SingularPrimes(R);
    wks:=[];
    for P in pp do
        _,p,_:=IsPrimePower(Index(R,P));
        k:=Valuation(Index(O,R),p);
        RP:=Order(ZBasis(R) cat ZBasis(O!!P^k));
        time Append(~wks,#WKICM_OLD(RP));
    end for;
    assert &*wks eq 114492;

*/
