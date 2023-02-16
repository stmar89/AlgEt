/* vim: set syntax=magma :*/

//freeze;

declare verbose OverOrders,3;

/////////////////////////////////////////////////////
// OverOrders for Etale Q algebras
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
// We are grateful to Edgar Costa (MIT) who helped 
// with a preliminary version of this code
/////////////////////////////////////////////////////

import "Ord.m" : crQZ , crZQ , Columns , hnf , MatrixAtoQ , MatrixAtoZ , MatrixQtoA , meet_zbasis , inclusion_matrix;

declare attributes AlgEtQOrd : MinimalOverOrders, // a sequence of tuples <P,{@ T1,...,Tn @}>, where
                                                  // P is a singular prime of the order R, and
                                                  // T1,..,Tn are all the minimal overorders with conductor (R:Ti)=P.
                               OverOrdersAtPrimes, // a sequence of tuples <P,{@ T1,...,Tn @}>, where
                                                  // P is a singular prime of the order R, and
                                                  // T1,..,Tn are all the overorders with P-primary conductor (R:Ti), plus R
                               OverOrders;  // all overorders


intrinsic IsMaximalAtPrime(R::AlgEtQOrd, P::AlgEtQIdl) -> BoolElt
{Returns whether R is maximal at the prime P, that is, if (R:O) is not contained in P, where O is the maximal order.}
    if IsMaximal(R) then
        return true;
    end if;
    if assigned R`SingularPrimes then
        return not P in SingularPrimes(R);
    end if;
    return not (Conductor(R) subset P);
end intrinsic;

intrinsic MinimalOverOrdersAtPrime(R::AlgEtQOrd, P::AlgEtQIdl) -> SetIndx[AlgEtQOrd]
{Given an order R and prime P of R, it returns the minimal overorders S of R with conductor (R:S) which is P-primary. The minimality assumption forces the conductor (R:S) to be exactly P. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana.}
    require IsPrime(P) : "The ideal P must be prime.";
    require Order(P) eq R : "P must be an ideal of R";

    if assigned R`MinimalOverOrders and exists(output){ tup : tup in R`MinimalOverOrders | tup[1] eq P } then
    // early exit if already computed
        output := output[2];
        return output;
    end if;

    // not already_computed
    output:={@ Universe({@ R @}) | @}; //empty set
    if not IsMaximalAtPrime(R,P) then
        zbR := ZBasis(R);
        // A minimal P-overorder S of R is contained in T=(P:P) and S/P is either 
        // a 2-dimensional vector space over R/P,
        // or a field extension of prime degree over R/P.
        // We call the overorder S of R contained in T that satisfy one of these two 
        // conditions 'potential minimal P-overoerders'i of R.
        pot_min_oo := {@ @};   // will contain all potential minimal P-overorders.
        pot_min_oo_2 := {@ @}; // will contain all potential minimal P-overorders where 
                               // the second necessary condition is satisfied
        F, f := ResidueField(P);      // f:R->R/P=F
        T := MultiplicatorRing(P);    // T=(P:P)
        V,mTV := QuotientVS(T, R, P); // mTV: T->T/R=V
        assert2 forall{ v : v in Basis(V) | mTV((v@@mTV)^2) in V };
        assert2 forall{ v : v in Basis(V) | mTV((v@@mTV)) eq v };
        assert2 forall{ t : t in ZBasis(T) | t-((mTV(t))@@mTV) in R };
        d := Dimension(V);
        if d eq 1 then // If d=1 then T/P had dimension 2 over R/T. This means that T is a minimal P-overorders of R
                       // See Proposition 5.3 of referenced paper
            Include(~output, T);
        else
            // minimal P-overorder S of R such that S/R has dimension 1 over R/P satisfy the followin property:
            // S/R is necessarily contained in an eigenspace of x->x^q acting on V=T/R.
            q:=#F;
            qpow:=hom<V->V | [mTV((v@@mTV)^q) : v in Basis(V)]>;
            eigen_vals:=[e[1] : e in Setseq(Eigenvalues(Matrix(qpow)))];
            eigen_spaces:=[Kernel(hom<V->V | [qpow(v)-e*v : v in Basis(V)]>): e in eigen_vals]; // eigenspaces in V
            subs_1:=[ W: W in &cat[Submodules(E) : E in eigen_spaces] | Dimension(W) eq 1];
            for W in subs_1 do
                // for each W of dim 1 we check whether is an order, that is, multiplicatively closed
                wT:=W.1@@mTV;
                if q eq 2 or mTV(wT^2) in W then
                    // for p eq 2 being a subspace of the eigenspace garantees that it is mult closed
                    S:=Order([wT] cat zbR : CheckIsKnownOrder:=false );
                    Include(~pot_min_oo,S);
                    Include(~output,S);// necessarly minimal, because it has dim 1
                end if;
            end for;
            // the other minimal overorders S of R are such that S/P is a finite field extension of prime degree of R/P
            dims := PrimesUpTo(d+1); //the plus one is to prevent issues when d=2.
            subs_2 := Submodules(V : CodimensionLimit := d-2); //we exclude dim 0 and 1
            subs_2 := [W : W in subs_2 | Dimension(W)+1 in dims]; // the +1 is there because we are 
                                                                  // working in T/R instead of T/P
            for W in subs_2 do //dim at least 2
                S := Order([(w@@mTV) : w in Basis(W)] cat zbR : CheckIsKnownOrder:=false );
                Include(~pot_min_oo,S);
                Include(~pot_min_oo_2,S);
            end for;
            //we remove non-minimals overorders from the potential ones in pot_min_oo_2
            for S in pot_min_oo_2 do
                if not exists {T : T in pot_min_oo | S ne T and T subset S} then
                    Include(~output, S);
                end if;
            end for;
        end if;
        // we check if any of these orders was already computed
        assert2 forall{S : S in output | ColonIdeal(R,R!!OneIdeal(S)) eq P};
        for i in [1..#output] do
            S:=output[i];
            IsKnownOrder(~S);
            ZBasisLLL(S);
        end for;
        if not assigned R`MinimalOverOrders then
            R`MinimalOverOrders:=[<P,output>];
        else
            Append(~R`MinimalOverOrders,<P,output>);
        end if;
    end if;
    return output;
end intrinsic;


intrinsic MinimalOverOrders(R::AlgEtQOrd) -> SetIndx[AlgEtQOrd]
{Computes the minimal overorders of R.}
    output:={@ Universe({@ R @}) | @}; //empty set
    if IsMaximal(R) then
        return output;
    end if;
    // Note: every overorders is a P-MinimalOverOrder for some singular prime P.
    pp:={@ P : P in SingularPrimes(R) @};
    if assigned R`MinimalOverOrders then
        done:={@ tup[1] : tup in R`MinimalOverOrders @};
        pp:=pp diff done;
    end if;
    for P in pp do 
        _:=MinimalOverOrdersAtPrime(R,P); //this populates the attribute R`MinimalOverOrders
    end for;
    output join:=&join[ tup[2] : tup in R`MinimalOverOrders ];
    return output; 
end intrinsic;

intrinsic OverOrdersAtPrime(R::AlgEtQOrd, P::AlgEtQIdl) -> SeqEnum[AlgEtQOrd]
{Given an order R and prime P of R, it returns R and the overorders S of R with conductor (R:S) which is P-primary. We recursively produce the minimal PP-overorders where PP are primes above P. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana.}
    require IsPrime(P) : "The ideal P must be prime.";
    require Order(P) eq R : "P must be an ideal of R";

    if assigned R`OverOrdersAtPrimes and exists(output){ tup : tup in R`OverOrdersAtPrimes | tup[1] eq P } then
    // early exit if already computed
        output := output[2];
        return output;
    end if;

    if assigned R`OverOrders then
    // early exit if already computed. This is useful when Loading the data using LoadWKICM
        return [ R ] cat [ S : S in R`OverOrders | PrimesAbove(ColonIdeal(R,R!!OneIdeal(S))) eq [ P ] ];
    end if;

    if IsMaximalAtPrime(R,P) then
        return [ R ];
    end if;

    ppO:=PrimesAbove(MaximalOrder(Algebra(R))!!P);
    queue := {@ R @};
    output:={@ R @};
    done:={@ @};
    while #queue gt 0 do
        pot_new:={@ @};
        for T in queue do
            pp:={@ OneIdeal(T) meet T!!Q : Q in ppO @};
            for i in [1..#pp] do
                Q:=pp[i];
                Q`IsPrime:=true;
            end for;
            //pp:=PrimesAbove(T!!P);
            for Q in pp do
                pot_new join:=MinimalOverOrdersAtPrime(T,Q);
            end for;
        end for;
        output join:=pot_new;
        done join:=queue;
        queue := pot_new diff done;
    end while;
    for iS in [1..#output] do
        S:=output[iS];
        ZBasisLLL(S);
    end for;
    output:=Setseq(output);
    assert2 forall{S : S in output | S eq R or PrimesAbove(ColonIdeal(R,R!!OneIdeal(S))) eq [ P ]};
    if not assigned R`OverOrdersAtPrimes then
        R`OverOrdersAtPrimes:=[<P,output>];
    else
        Append(~R`OverOrdersAtPrimes,<P,output>);
    end if;
    return output;
end intrinsic;

intrinsic OverOrders(R::AlgEtQOrd : populateoo_in_oo:=false) -> SeqEnum[AlgEtQOrd]
{We compute all the overorders of R. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana. The Vararg "populateoo_inoo" (default false) determines whether we should fill the attribute T`OverOrders for every overorder T of R.}
    if not assigned R`OverOrders then
        if IsMaximal(R) then
            output:=[R];
        else
            output:=[];
            pp:=SingularPrimes(R);
            vtime OverOrders,2 : oo_at_Ps:=[ OverOrdersAtPrime(R,P) : P in pp ];
            vprintf OverOrders,2 : "P-overorders for each P %o\n",[#x : x in oo_at_Ps];
            vprintf OverOrders,2 : "total number of overorders %o\n",&*[#x : x in oo_at_Ps];
            assert forall{x : x in oo_at_Ps | x[1] eq R};
            cc:=CartesianProduct([[1..#x] : x in oo_at_Ps]);
            for c in cc do
                S:=[ c[j] : j in [1..#c] ];
                // R is contained in all the orders. 
                // Using this info we can make the generation process a bit faster.
                if forall{x : x in S | x eq 1 } then
                    S:=R;
                elif #[x : x in S | x ne 1] eq 1 then
                    // this case happens for the minimal overorders of R, 
                    // and when there is only one singular prime of R
                    _,j:=Max(S); // the position of the non-1 entry
                    S:=oo_at_Ps[j][S[j]]; // the only order in the tuple bigger than R
                else
                    gens:=&cat[ZBasis(oo_at_Ps[j][S[j]]) : j in [1..#S] | S[j] ne 1];
                    S:=Order(gens);
                end if;
                ZBasisLLL(S);
                Append(~output,S);
            end for;
            assert2 R in output;
            assert2 MaximalOrder(Algebra(R)) in output;
        end if;
        R`OverOrders:=output;
    end if;
    // there might be a better way to do this
    if populateoo_in_oo then
        for i in [1..#R`OverOrders] do
            S := R`OverOrders[i];
            if not assigned S`OverOrders then
                S`OverOrders := [ T : T in R`OverOrders | S subset T ];
            end if;
        end for;
    end if;
    return R`OverOrders;
end intrinsic;

intrinsic FindOverOrders(R::AlgEtQOrd : populateoo_in_oo:=false) -> SetIndx[AlgEtQOrd]
{We compute all the overorders of R. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana. The Vararg "populateoo_inoo" (default false) determines whether we should fill the attribute T`OverOrders for every overorder T of R.}
    return OverOrders(R : populateoo_in_oo:=populateoo_in_oo);
end intrinsic;

/* TESTS

    printf "### Testing OverOrders:";
	AttachSpec("~/packages_github/AlgEt/spec");

    SetVerbose("OverOrders",1);
    SetAssertions(1);

    _<x>:=PolynomialRing(Integers());
    f:=(x^4+16);
    A:=EtaleAlgebra(f);
    O:=MaximalOrder(A);
    assert #FindOverOrders(O) eq 1;
    assert #MinimalOverOrders(O) eq 0;

    E:=EquationOrder(A);
    oo:=FindOverOrders(E);
    assert #oo eq 11;

    printf " all good!\n"; 

*/
