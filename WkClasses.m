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

intrinsic WKICM_bar(S::AlgEtOrd) -> SeqEnum
{returns all the weak eq classes I, such that (I:I)=S}
    if not assigned S`WKICM_bar then
        if IsGorenstein(S) then
            S`WKICM_bar:=[OneIdeal(S)];
        else
            if CohenMacaulayType(S) eq 2 then
                S`WKICM_bar:=wkicm_bar_CM_type2(S,NonGorensteinPrimes(S));
            else
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
                //T_ZBasis:=ZBasis(T);
                ff:=ColonIdeal(S,S!!OneIdeal(T)); //the relative conductor (S:T)
                gens_ff_over_S:=Generators(ff);
                ffT:=T!!ff; //ff is a T module
                //OLD CODE
                //ff_ZBasis:=ZBasis(ff);
                //seqWk_bar:=[];
                //F:=FreeAbelianGroup(Dimension(Algebra(S)));
                //matT:=Matrix(T_ZBasis);
                //matff:=Matrix(ff_ZBasis);
                //rel:=[F ! Eltseq(x) : x in Rows(matff*matT^-1)];
                //Q,q:=quo<F|rel>; //Q=T/(S:T)
                Q,q:=ResidueRing(T,ffT); //Q=T/(S:T) , q:T->Q
                vprintf WkClasses,2:"T/(S:T) = %o\n",Q;
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
                S`WKICM_bar:=seqWk_bar;
            end if;
        end if;
    end if;
    return S`WKICM_bar;
end intrinsic;

intrinsic WKICM(E::AlgEtOrd)->SeqEnum
{computes the Weak equivalence class monoid of E}
    if not assigned E`WKICM then
        seqOO:=FindOverOrders(E : populateoo_in_oo:=true);
        E`WKICM:=&cat[[(E!!I) : I in WKICM_bar(S)] : S in seqOO ];
    end if;
    return E`WKICM;
end intrinsic;

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




