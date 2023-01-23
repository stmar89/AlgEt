/* vim: set syntax=magma :*/

freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose new_wk_icm_bar, 2;

//------------
// NOTES: in the end we will need ideals up to isomorphism. I think that in the recursion we can skip the 1-dimensional vector spaces (which are probably the wast majority) since they will all give rise to the invertible weak equivalence class. Need to think about it more.
//------------

// when incorporating in package, the next line will not be necessary
import "../AlgEtQ/LowCohenMacaulayType.m" : wkicm_bar_CM_type2;

intrinsic IntermediateIdealsVSWithTrivialExtensionAndPrescribedMultiplicatorRing(I::AlgEtQIdl,P::AlgEtQIdl, O::AlgEtQOrd)->SetIndx[AlgEtQIdl]
{Given a fractional S-ideal I, a prime ideal P of S, and an overorder O of S such that O subset (I:I), it returns all the fractional S-ideals J such that PI subset J subset I, (J:J) eq S, and JO=I.}
    require IsPrime(P) : "The ideal P must be prime";
    S:=Order(P);
    require S eq Order(I) : "The ideals must be over the same order";
    require S subset O : "O is not an overorder of Order(I)";
    IO:=O!!I;
    require O subset MultiplicatorRing(I) : "I is not an O-ideal";
    PI:=P*I;
    zbPI:=ZBasis(PI);
    Q,q:=QuotientVS(I,P*I,P); // q:I->I/PI=Q
    zbO:=ZBasis(O);

    queue:={@ Q @};
    output:={@  @};
    if MultiplicatorRing(I) eq S then
        Include(~output, S!!I);
    end if;
    done:={@ @};
    while #queue gt 0 do
        pot_new:={@ x : x in &cat[ MaximalSubmodules(elt) : elt in queue ]@};
        pot_new_lifts:=[ Ideal(S,[ (Q!b)@@q : b in Basis(W) ] cat zbPI) : W in pot_new]; 
        indices:=[ i : i in [1..#pot_new] | O!!pot_new_lifts[i] eq IO ]; 
            // checking trivial ext with the lifts it is faster than doing it in Q
        pot_new:={@ pot_new[i] : i in [1..#pot_new] | i in indices @};
            // only the ones with trivial extension
        pot_new_lifts:={@ pot_new_lifts[i] : i in [1..#pot_new_lifts] | i in indices @};
        output join:={@ pot_new_lifts[i] : i in [1..#pot_new] | MultiplicatorRing(pot_new_lifts[i]) eq S @};
            // one the lifts with multiplicator ring S
        done join:=queue;
        // If W does not have trivial extension, then all its subspace will not have trivial extension.
        // Hence we do not add such W to the queue.
        queue := pot_new diff done; 
    end while;
    return output;
end intrinsic;

intrinsic WKICM_bar(S::AlgEtQOrd : Method:="Auto") -> SeqEnum
{Returns all the weak eq classes I, such that (I:I)=S. The VarArg Method (default "Auto") determines if we should use the "IntermediateIdeals" routine or the "LowIndexProcess", which is potentially much slower but more memory efficient.}
    if not assigned S`WKICM_bar then
        if IsGorenstein(S) then
            vprintf new_wk_icm_bar,2 : " Gorenstein case\n";
            S`WKICM_bar:=[OneIdeal(S)];
        elif CohenMacaulayType(S) eq 2 then
            vprintf new_wk_icm_bar,2 : " Cohen Macaulay type 2 case\n";
            S`WKICM_bar:=wkicm_bar_CM_type2(S,NonGorensteinPrimes(S));
        else
            require Method in {"Auto","LowIndexProcess","IntermediateIdeals","IntermediateIdealsVSWithTrivialExtensionAndPrescribedMultiplicatorRing"} : "The VarArg parameter Method is assigned to a not avaialble value";
            vprintf WkClasses,2:"Order of CohenMacaulayType = %o\n",CohenMacaulayType(S);
            // general case
            seqWk_bar:=[];
            if Method in {"Auto","IntermediateIdealsVSWithTrivialExtensionAndPrescribedMultiplicatorRing"} then
                vprint new_wk_icm_bar,2 : "Using new method";
                pp:=PrimesAbove(Conductor(S));
                oo:=FindOverOrders(S);
                mult_pp:=[ oo[Index(oo,MultiplicatorRing(P))] : P in pp ];
                assert forall{T : T in mult_pp | assigned T`WKICM_bar};
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
                P:=pp[iP];
                T:=mult_pp[iP];
                wkT:=WKICM_bar(mult_pp[iP]);
                cands:=&join[ IntermediateIdealsVSWithTrivialExtensionAndPrescribedMultiplicatorRing(S!!J,P,T) : J in wkT ];
                for I in cands do
                  if not exists{J : J in seqWk_bar | IsWeakEquivalent(I,J)} then
                    ZBasisLLL(I);
                    Append(~seqWk_bar,I);
                  end if;
                end for;
            else
              St:=TraceDualIdeal(S);
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
    "NEW method";
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
    require Method in {"Auto","LowIndexProcess","IntermediateIdeals","IntermediateIdealsVSWithTrivialExtensionAndPrescribedMultiplicatorRing"} : "The VarArg parameter Method is assigned to a not available value";
    if not assigned E`WKICM then
        wk:=[];
        oo:=FindOverOrders(E : populateoo_in_oo:=true);
        if Method in {"Auto","IntermediateIdealsVSWithTrivialExtensionAndPrescribedMultiplicatorRing"} then
            // We want to compute WKICM_bar(S) for S running from top to bottom 
            // of the lattice of overorder of R.
            // The reason is that when we want to compute it for S, if S is not Gorenstein or of CMType 2, then
            // we need to have WKICM_bar(T) computed for T=(P:P) for all P singular prime of S.
            // Hence we proceed recursively starting from O, and then filling it for the Maximal SubOrders.
            O:=MaximalOrder(Algebra(E));
            edges:=Edges(GraphOverOrders(E));
            edges:=[ [TerminalVertex(e),InitialVertex(e)] : e in edges ]; // we reverse all the inclusion
            // D is the directed graph of (reverse) minimal inclusion of overorders of R.
            D:=Digraph<#oo|edges>;
            vprint new_wk_icm_bar,2 : "We have computed the directed graph of (reverse) inclusion of overorders of R.";

            vv:=VertexSet(D);
            max_dist:=Distance(vv.Index(oo,O),vv.Index(oo,E));
            for d in [0..max_dist] do
                vprintf new_wk_icm_bar,2 : "Dealing now with orders of distance %o from the maximal order.\n",d;
                oo_d:=Sphere(vv.Index(oo,O),d);
                for v in oo_d do
                    vprintf new_wk_icm_bar,2 : "Computing WKICM_bar for the %o-th order ... \n",v;
                    T:=oo[Index(v)];
                    vtime new_wk_icm_bar,2 :wk cat:=[(E!!I) : I in WKICM_bar(T : Method:=Method)];
                end for;
            end for;
        else
            wk:=&cat[[(E!!I) : I in WKICM_bar(S : Method:=Method)] : S in oo ];
        end if;
        for I in wk do
            ZBasisLLL(I);
        end for;
        E`WKICM:=wk;
    end if;
    return E`WKICM;
end intrinsic;


/* TESTS
 
    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/dev/new_wk_icm.m");
    P<x>:=PolynomialRing(Integers());
    f:=x^8+16;
    A:=EtaleAlgebra(f);
    R:=EquationOrder(A);
    time _:=FindOverOrders(R); // <4 secs
    time assert #WKICM(R) eq 173; // ~40 secs
    "the size of the output, 173 classes, has been computed using the OLD method below, in 189000 seconds";
    // quit;

//    magma
// 
//     AttachSpec("~/packages_github/AlgEt/spec");
//     "OLD method";
//   	_<x>:=PolynomialRing(Integers());
//     f:=x^8+16; 
//     A:=EtaleAlgebra(f);
//     R:=EquationOrder(A);
//     time _:=FindOverOrders(R);
//     time #WKICM(R);
//     quit;
//     screen -r test_new_wk_icm
// 
//     quit;
//     cd ~/packages_github/AlgEt/
//     git pull; sleep 1;
//     magma


    // looping over slow input
    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/dev/new_wk_icm.m");
    P<x>:=PolynomialRing(Integers());
    str:=Split(Read("~/packages_github/AlgEt/dev/input_slow_sorted.txt"));
    poly:=[ P!eval(l) : l in str[1..10] ];
    for i->f in poly do
        A:=EtaleAlgebra(f);
        F:=PrimitiveElement(A);
        q:=Integers() ! Round(ConstantCoefficient(f)^(2/Degree(f)));
        R:=Order([F,q/F]);
        _:=FindOverOrders(R);
        t0:=Cputime();
        printf "i=%o f=%o #wk=%o seconds=%o\n",i,f,#WKICM(R),Cputime(t0);
    end for;

    // freakking big
    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/dev/new_wk_icm.m");
    SetVerbose("new_wk_icm_bar",2);
    P<x>:=PolynomialRing(Integers());
    f:=x^8 - 2*x^7 + 7*x^6 - 14*x^5 + 40*x^4 - 56*x^3 + 112*x^2 - 128*x + 256;
    A:=EtaleAlgebra(f);
    R:=EquationOrder(A);
    "Computing overorders ...";
    time oo:=FindOverOrders(R); //14824 seconds, 16320 overorders
    #oo, "overorders found.";
    wk:=WKICM(R); // 114492

    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/dev/new_wk_icm.m");
    //SetVerbose("new_wk_icm_bar",2);
    P<x>:=PolynomialRing(Integers());
    f:=x^8 - 2*x^7 + 7*x^6 - 14*x^5 + 40*x^4 - 56*x^3 + 112*x^2 - 128*x + 256;
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
    assert &*wks eq 114492;

// ISSUE with the next: Segmentation fault for the Equation Order. Maybe there are just too many subscpaces
    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/dev/new_wk_icm.m");
    SetVerbose("new_wk_icm_bar",2);
    P<x>:=PolynomialRing(Integers());
    f:=x^6 + 8*x^5 + 50*x^4 + 200*x^3 + 1250*x^2 + 5000*x + 15625;
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    // R:=Order([F,25/F]);
    R:=EquationOrder(A); //3200 overorders
    "Computing overorders ...";
    time oo:=FindOverOrders(R);
    #oo, "overorders found.";
    wk:=WKICM(R);

*/
