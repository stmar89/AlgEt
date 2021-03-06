/* vim: set syntax=magma :*/

//freeze;

/////////////////////////////////////////////////////
// Ideal class monoid and weak equivalence classes for orders in Etale Q algebras
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose WkClasses, 3;

declare attributes AlgEtOrd:WKICM,
                            WKICM_bar;

import "LowCohenMacaulayType.m" : wkicm_bar_CM_type2;

intrinsic WKICM_bar(S::AlgEtOrd : Method:="Auto") -> SeqEnum
{ Returns all the weak eq classes I, such that (I:I)=S. The VarArg Method (default "Auto") determines if we should use the "IntermediateIdeals" routine or the "LowIndexProcess", which is potentially much slower but more memory efficient.}
    if not assigned S`WKICM_bar then
        if IsGorenstein(S) then
            S`WKICM_bar:=[OneIdeal(S)];
        else
            if CohenMacaulayType(S) eq 2 then
                S`WKICM_bar:=wkicm_bar_CM_type2(S,NonGorensteinPrimes(S));
            else
                require Method in {"Auto","LowIndexProcess","IntermediateIdeals"} : "The VarArg parameter Method is assigned to a not avaialble value";
                vprintf WkClasses,2:"Order of CohenMacaulayType = %o\n",CohenMacaulayType(S);
                // general case
                //TODO : prime per prime;
                A:=Algebra(S);
                degA:=Dimension(A);
                seqWk_bar:=[];
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
                     
                    // NEW CODE based on IntermediateIdealsWithPrescribedMultiplicatorRing
                    vprintf WkClasses,2:"Using code with IntermediateIdealsWithPrescribedMultiplicatorRing: T/(S:T) = %o\n",Q;
                    idls:=IntermediateIdealsWithPrescribedMultiplicatorRing(T1,ff);
                    for I in idls do
                        if not exists{ J : J in seqWk_bar | IsWeakEquivalent(I,J) } then
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
                        //geninF:=[(f(QP ! x))@@q : x in Generators(H)];
                        //coeff:=[Eltseq(x) : x in geninF];
                        //I:=ideal<S| [&+[T_ZBasis[i]*x[i] : i in [1..degA]] : x in coeff] cat ff_ZBasis>;
                        I:=Ideal(S, [ (f(QP!x))@@q : x in Generators(H) ] cat gens_ff_over_S);
                        if not I in seqWk_bar and 
                            MultiplicatorRing(I) eq S and 
                            not exists{J : J in seqWk_bar | IsWeakEquivalent(I,J)} then 
                                Append(~seqWk_bar,I);
                        end if;
                    end while;
                end if;
                S`WKICM_bar:=seqWk_bar;
            end if;
        end if;
    end if;
    return S`WKICM_bar;
end intrinsic;

intrinsic WKICM(E::AlgEtOrd : Method:="Auto")->SeqEnum
{Computes the Weak equivalence class monoid of E. The VarArg Method (default "Auto") determines if we should use the "IntermediateIdeals" routine or the "LowIndexProcess", which is potentially much slower but more memory efficient.}
    if not assigned E`WKICM then
        require Method in {"Auto","LowIndexProcess","IntermediateIdeals"} : "The VarArg parameter Method is assigned to a not avaialble value";
        seqOO:=FindOverOrders(E : populateoo_in_oo:=true);
        E`WKICM:=&cat[[(E!!I) : I in WKICM_bar(S : Method:=Method)] : S in seqOO ];
    end if;
    return E`WKICM;
end intrinsic;

/* 
intrinsic WKICM_bar_intermediate_idls(S::AlgEtOrd) -> SeqEnum[AlgEtIdl]
{???}
// add attributes ???
//  to be used only for non Gorenstein orders
    St:=TraceDualIdeal(S);
    T:=&meet([ T : T in FindOverOrders(S) | IsInvertible(T !! St) ]);
    //this construction of T is conjectural, hence the next assert. If the assert fails, please report it.
    assert IsInvertible(T !! St);
    return output;
end intrinsic;
*/

/*TEST

	AttachSpec("~/packages_github/AlgEt/spec");

	SetAssertions(2);

	_<x>:=PolynomialRing(Integers());
    f:=x^4-1000*x^3-1000*x^2-1000*x-1000;
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    SeqWC:=WKICM(E);
    if #SeqWC ne 25 then
      test:=false;
      printf"\nERROR: SeqWC of f=%o\n",f;
    end if;
    // _:=ICM(E); //computing the Pics is very slow!

    f:=x^4+291*x^3-988*x^2-1000*x-1000;
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    SeqWC:=WKICM(E);
    if #SeqWC ne 20 then
      test:=false;
      printf"\nERROR: SeqWC of f=%o\n",f;
    end if;
    _:=ICM(E);

    f:=x^3+31*x^2+43*x+77;
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    if #FindOverOrders(E) ne 15 then 
      test:=false;
      printf"\nERROR: OverOrders of f=%o\n",f;
    end if;
    SeqWC:=WKICM(E);
    if #SeqWC ne 23 then
      test:=false;
      printf"\nERROR: SeqWC of f=%o\n",f;
    end if;
    _:=ICM(E);

*/




