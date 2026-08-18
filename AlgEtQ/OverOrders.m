/////////////////////////////////////////////////////
// Copyright 2025.
// Stefano Marseglia, stefano.marseglia89@gmail.com
// https://stmar89.github.io/index.html
//
// Distributed under the terms of the CC-BY 4.0 licence.
// https://creativecommons.org/licenses/by/4.0/
/////////////////////////////////////////////////////


// Copyright 2023, S. Marseglia
/////////////////////////////////////////////////////

freeze;

declare verbose OverOrders, 1;

declare attributes AlgEtQOrd : MinimalOverOrders, // a sequence of tuples <P,{@ T1,...,Tn @}>, where
                                                  // P is a singular prime of the order R, and
                                                  // T1,..,Tn are all the minimal overorders with conductor (R:Ti)=P.
                               OverOrdersAtPrimes, // a sequence of tuples <P,{@ T1,...,Tn @}>, where
                                                  // P is a singular prime of the order R, and
                                                  // T1,..,Tn are all the overorders with P-primary conductor (R:Ti), plus R
                               OverOrders;  // all overorders

///# Overorders
/// Let $R$ be an order in an étale algebra $A$ over $\mathbb{Q}$ with maximal order $\mathcal{O}_A$.
/// Since the quotient $\mathcal{O}_A/R$ is finite, it follows that there are only finitely many overorders of $R$.
/// The intrinsic `OverOrders` is an implementation of the Hofmann-Sircana algorithm contained in "On the computation of overorders" (Int. J. Number Theory 16, No. 4, 857-879 (2020)).
/// Recall that if $\mathfrak{p}$ is a prime of $R$, then an overorder $S$ of $R$ is a $\mathfrak{p}$-overorder of $R$ if the colon ideal $(R:S)$ is a $\mathfrak{p}$-primary ideal, or, equivalently, the finite $R$-module $S/R$ is only supported at $\mathfrak{p}$.
/// The algorithm can be summarized in two steps.
/// The first one builds on the observation that the lattice of inclusions of the overorders of $R$ is the product of the lattice of inclusions of the $\mathfrak{p}$-overorders of $R$, where $\mathfrak{p}$ runs over the singular primes of $R$.
/// The second step consists of computing the $\mathfrak{p}$-overorders of $R$ by recursively constructing minimal $\mathfrak{p}$-overorders.

/// Returns whether the order $R$ is maximal at the prime $P$, that is, if $(R:O)$ is not contained in $P$, where $O$ is the maximal order.
intrinsic IsMaximalAtPrime(R::AlgEtQOrd, P::AlgEtQIdl) -> BoolElt
{Returns whether the order R is maximal at the prime P, that is, if (R:O) is not contained in P, where O is the maximal order.}
    if IsMaximal(R) then
        return true;
    end if;
    if assigned R`SingularPrimes then
        return not P in SingularPrimes(R);
    end if;
    return not (Conductor(R) subset P);
end intrinsic;

/// Given an order $R$ and prime $P$ of $R$, returns the minimal overorders $S$ of R whose conductor $(R:S)$ is a $P$-primary ideal of $R$. The minimality assumption forces the conductor $(R:S)$ to be exactly $P$.
intrinsic MinimalOverOrdersAtPrime(R::AlgEtQOrd, P::AlgEtQIdl) -> SeqEnum[AlgEtQOrd]
{Given an order R and prime P of R, returns the minimal overorders S of R whose conductor (R:S) is a P-primary ideal of R. The minimality assumption forces the conductor (R:S) to be exactly P. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana.}
    require IsPrime(P) : "The ideal P must be prime.";
    require Order(P) eq R : "P must be an ideal of R";

    if assigned R`MinimalOverOrders and exists(output){ tup : tup in R`MinimalOverOrders | tup[1] eq P } then
    // early exit if already computed
        output := output[2];
        return output;
    end if;

    output:=[ Universe([ R ]) | ]; //empty sequence
    if not IsMaximalAtPrime(R,P) then
        zbR := ZBasis(R);
        // A minimal P-overorder S of R is contained in T=(P:P), and S/P is either
        // 1) a 2-dimensional vector space over R/P,
        // 2) or a field extension of prime degree >2 over R/P.
        // These two conditions are mutually exclusive.
        min_oo_type1 := [ ];   // will contain minimal P-overorders of the first type
        min_oo_type2 := [ ]; // will contain all potential minimal P-overorders where
                               // the second necessary condition is satisfied
        T := MultiplicatorRing(P);     // T=(P:P)
        VT,mTtoVT:= QuotientVS(T,P,P); // VT=T/P is a R/P-algebra.
        dVT := [1..Dimension(VT)];
        if #dVT eq 2 then
            // VT is 2-dimensional if and only if T is the unique minimal P-overorder.
            Append(~output, T);
        else
            VT1 := mTtoVT(1);
            V,mVTtoV:= quo<VT|[VT1]>;      // V is not a ring since it does not have a well defined multiplication.
            dV := Dimension(V);
            assert #dVT eq dV+1; 
            // minimal P-overorder S of R such that S/R has dimension 1 over R/P satisfy the following property:
            // S/R is necessarily contained in an eigenspace of x->x^q acting on V=T/R.
            mTV := map<T->V|x:->mVTtoV(mTtoVT(x)),y:->y@@mVTtoV@@mTtoVT>; // mTV: T->T/R=V
            // We create 3 auxiliary function: a multiplication map on VT, a n-power map on VT, and 
            // a functions that checks whether a subspace of V is mult closed, that is, the reduction of an order. 
            // The latter builds upon the first 2.
            // This approach keeps the coefficients bounded, allowing to use linear algebra over finite fields 
            // rather than over the integers.
            VTi_inT := [VT.i@@mTtoVT:i in dVT];
            VT0:=Zero(VT);
            VTi_VTj := [[VT0:i in dVT]:j in dVT];
            for j in dVT,i in [1..j] do
                x:=mTtoVT(VTi_inT[i]*VTi_inT[j]);
                VTi_VTj[i,j]:=x;
                VTi_VTj[j,i]:=x;
            end for;
            mult := func<x,y|&+[x[i]*y[j]*VTi_VTj[i,j]:i,j in dVT]>; //x*y in VT
            pow:=function(x,n)
                // x^n in VT
                result := VT1;
                base := x;
                exp := n;
                while exp gt 0 do
                    if exp mod 2 eq 1 then
                        result := mult(result,base);
                    end if;
                    base := mult(base,base);
                    exp := exp div 2;
                end while;
                return result;
            end function;
            is_ord_inV:=function(W)
                // Given a subvectors space W of V, returns whether its preimage in VT is multiplicatively closed.
                // This happens precisely when W is the reduction of an order S such that R < S < T.
                dWT:=Dimension(W)+1;
                bWTinVT:=Append([g@@mVTtoV:g in Basis(W)],VT1);
                return forall{i:i in [1..dWT]|forall{j:j in [1..i]|mVTtoV(mult(bWTinVT[i],bWTinVT[j])) in W}};
            end function;
            assert2 forall{z:z in zbR|mTV(z) eq 0};
            assert2 forall{ v : v,w in Basis(VT) | mult(v,w) eq mTtoVT(v@@mTtoVT * w@@mTtoVT) };
            assert2 forall{ v : v in Basis(VT) | pow(v,10) eq mTtoVT((v@@mTtoVT)^10)};
            assert2 forall{ v : v in Basis(VT) | mTV((v@@mTtoVT)^2) in VT };
            assert2 forall{ v : v in Basis(VT) | pow(v,2) in VT };
            assert2 forall{ v : v,w in Basis(VT) | mult(v,w) in VT };
            assert2 forall{ v : v in Basis(VT) | mTtoVT((v@@mTtoVT)) eq v };
            assert2 forall{ t : t in ZBasis(T) | t-((mTtoVT(t))@@mTtoVT) in R };

            q:=Index(R,P);
            qpow:=hom<V->V | [mVTtoV(pow(v@@mVTtoV,q)) : v in Basis(V)]>;
            eigen_vals:=[e[1] : e in Setseq(Eigenvalues(Matrix(qpow)))];
            eigen_spaces:=[Kernel(hom<V->V | [qpow(v)-e*v : v in Basis(V)]>): e in eigen_vals]; // eigenspaces in V

            //FIXME check which one of the next two lines is better/correct
            subs_1:=[ W: W in &cat[Submodules(E) : E in eigen_spaces] | Dimension(W) eq 1];
            //subs_1:=[ W: W in &cat[Submodules(E : CodimensionLimit:=Dimension(E)-1) : E in eigen_spaces] | Dimension(W) eq 1 ];
            for W in subs_1 do
                // for each W of dim 1 we check whether is an order, that is, multiplicatively closed
                if q eq 2 then
                    // for q eq 2 being a subspace of the eigenspace guarantees that it is mult closed
                    w:=V!W.1;
                    Append(~min_oo_type1,W);
                    S:=Order([(V!v)@@mTV:v in Basis(W)] cat zbR : CheckIsKnownOrder:=true );
                    Append(~output,S);// necessarily minimal, because it has dim 1
                else
                    wVT:=(V!W.1)@@mVTtoV;
                    if mVTtoV(pow(wVT,2)) in W then
                        Append(~min_oo_type1,W);
                        S:=Order(Append(zbR,wVT@@mTtoVT) : CheckIsKnownOrder:=true );
                        Append(~output,S);// necessarily minimal, because it has dim 1
                    end if;
                end if;
            end for;
            // the other minimal overorders S of R are such that S/P is a finite field extension of prime degree of R/P
            dims := PrimesUpTo(dV+1); //the plus one is to prevent issues when d=2.
            subs_type2 := Submodules(V : CodimensionLimit := dV-2); //we exclude dim 0 and 1
            for W in subs_type2 do
                if Dimension(W)+1 in dims and is_ord_inV(W) then
                    assert Dimension(W) ge 2;
                    S:=Order([v@@mTV:v in Basis(W)] cat zbR : CheckIsKnownOrder:=true );
                    Append(~output,S);
                end if;
            end for;
        end if;
        // we check if any of these orders was already computed
        assert2 forall{S : S in output | ColonIdeal(R,R!!OneIdeal(S)) eq P};
        if not assigned R`MinimalOverOrders then
            R`MinimalOverOrders:=[<P,output>];
        else
            Append(~R`MinimalOverOrders,<P,output>);
        end if;
    end if;
    return output;
end intrinsic;


/// Computes the minimal overorders of the order $R$.
intrinsic MinimalOverOrders(R::AlgEtQOrd) -> SeqEnum[AlgEtQOrd]
{Computes the minimal overorders of R.}
    output:=[ Universe([ R ]) | ]; //empty set
    if IsMaximal(R) then
        return output;
    end if;
    // Note: every minimal overorder is a P-MinimalOverOrder for a unique singular prime P.
    // i.e. if P ne Q, then the set of P-minimal overorders and Q-minimal overorders are disjoint!
    return &cat[ MinimalOverOrdersAtPrime(R,P):P in SingularPrimes(R) ];
// FIXME remove the following old code after testing
//    pp:={@ P : P in SingularPrimes(R) @};
//    if assigned R`MinimalOverOrders then
//        done:={@ tup[1] : tup in R`MinimalOverOrders @};
//        pp:=pp diff done;
//    end if;
//    for P in pp do
//        _:=MinimalOverOrdersAtPrime(R,P); //this populates the attribute R`MinimalOverOrders
//    end for;
//    output join:=&join[ tup[2] : tup in R`MinimalOverOrders ];
//    return output;
end intrinsic;


// Given an order $R$ and prime $P$ of $R$, returns $R$ and the overorders $S$ of $R$ with conductor $(R:S)$ which is $P$-primary.
intrinsic OverOrdersAtPrime(R::AlgEtQOrd, P::AlgEtQIdl) -> SeqEnum[AlgEtQOrd]
{Given an order R and prime P of R, returns R and the overorders S of R with conductor (R:S) which is P-primary. We recursively produce the minimal PP-overorders where PP are primes above P. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana.}
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
            for Q in pp do
                pot_new join:=SequenceToIndexedSet(MinimalOverOrdersAtPrime(T,Q));
            end for;
        end for;
        output join:=pot_new;
        done join:=queue;
        queue := pot_new diff done;
    end while;
    for iS in [1..#output] do
        S:=output[iS];
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

/// Returns the overorders of $R$. The parameter `populateoo_in_oo` (default false) determines whether the algorithms fills the attribute `OverOrders` for each overorder of $R$.
intrinsic OverOrders(R::AlgEtQOrd : populateoo_in_oo:=false) -> SeqEnum[AlgEtQOrd]
{Returns the overorders of R. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana. The parameter "populateoo_in_oo" (default false) determines whether we should fill the attribute T`OverOrders for every overorder T of R.}
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

///ditto
intrinsic FindOverOrders(R::AlgEtQOrd : populateoo_in_oo:=false) -> SetIndx[AlgEtQOrd]
{Returns the overorders of R. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana. The parameter "populateoo_in_oo" (default false) determines whether we should fill the attribute T`OverOrders for every overorder T of R.}
    return OverOrders(R : populateoo_in_oo:=populateoo_in_oo);
end intrinsic;

///# Example 5
/// ```
/// _<x>:=PolynomialRing(Integers());
/// f:=(x^4+16)*(x^4+81);
/// A:=EtaleAlgebra(f);
/// E:=EquationOrder(A);
/// oo:=OverOrders(E);
/// #oo;
/// pp:=SingularPrimes(E);
/// // We see that the size of the lattice of inclusions of the overorders is the product of the sizes of local component of the lattice.
/// #oo eq &*[ #OverOrdersAtPrime(E,P) : P in pp ];
/// // Now we consider only the P-overorders S for the first singular prime P.
/// // We verify that there is always a positive integer i such that (R:S)^i is invertible in its multiplicator ring.
/// ooP:=OverOrdersAtPrime(E,pp[1]);
/// #ooP;
/// for S in ooP do
///     C:=ColonIdeal(E,E!!OneIdeal(S));
///     [ IsInvertible(Ti!!Ci) where Ti:=MultiplicatorRing(Ci) where Ci:=C^i : i in [1..10]];
/// end for;
/// ```

/* TESTS

    printf "### Testing OverOrders:";
    //AttachSpec("~/packages_github/AlgEt/spec");

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

    printf " all good!";
*/
