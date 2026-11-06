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
                               OverOrdersAtPrimes, // a sequence of tuples <P,[ T1,...,Tn ]>, where
                                                  // P is a singular prime of the order R, and
                                                  // T1,..,Tn are all the overorders with P-primary conductor (R:Ti), plus R
                               OverOrders;        // all overorders

sort_overorders:=procedure(~oo)
// Given a sequence oo of orders, sorts it by decreasing Index(S), breaking ties using myHash(S).
// If the orders in oo are overorders of a fixed order R, then Index(S,R)=Index(R)/Index(S), so this is
// the same as sorting by increasing index [S:R]. In particular, R is always the first entry, and the
// maximal order the last one.
// This makes the output of OverOrdersAtPrime and OverOrders independent of the way the orders are produced.
    keys:=[ [1/Index(S)] cat myHash(S) : S in oo ];
    ParallelSort(~keys,~oo);
end procedure;

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

/// Given an order $R$ and prime $P$ of $R$, returns the minimal overorders $S$ of R whose conductor $(R:S)$ is a $P$-primary ideal of $R$. The minimality assumption forces the conductor $(R:S)$ to be exactly $P$. The vararg `CheckIsKnownOrder` (default `true`) determines whether the output is checked against the global database of known orders. If it is set to `false`, the output is not stored in the attribute `MinimalOverOrders` of $R$. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana.
intrinsic MinimalOverOrdersAtPrime(R::AlgEtQOrd, P::AlgEtQIdl : CheckIsKnownOrder:=true) -> SetIndx[AlgEtQOrd]
{Given an order R and prime P of R, returns the minimal overorders S of R whose conductor (R:S) is a P-primary ideal of R. The minimality assumption forces the conductor (R:S) to be exactly P. The vararg CheckIsKnownOrder (default true) determines whether the output is checked against the global database of known orders. If it is set to false, the output is not stored in the attribute MinimalOverOrders of R. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana.}
    require IsPrime(P) : "The ideal P must be prime.";
    require Order(P) eq R : "P must be an ideal of R";

    if assigned R`MinimalOverOrders and exists(output){ tup : tup in R`MinimalOverOrders | tup[1] eq P } then
    // early exit if already computed
        output := output[2];
        return output;
    end if;

    output:={@ Universe({@ R @}) | @}; //empty set
    if not IsMaximalAtPrime(R,P) then
        zbR := ZBasis(R);
        // A minimal P-overorder S of R is contained in T=(P:P), and S/P is either
        // 1) a 2-dimensional vector space over R/P,
        // 2) or a field extension of prime degree >2 over R/P.
        // These two conditions are mutually exclusive.
        T := MultiplicatorRing(P);     // T=(P:P)
        VT,mTtoVT:= QuotientVS(T,P,P); // VT=T/P is a R/P-algebra.
        dVT := [1..Dimension(VT)];
        if #dVT eq 2 then
            // VT is 2-dimensional if and only if T is the unique minimal P-overorder.
            Include(~output, T);
        else
            VT1 := mTtoVT(1);
            V,mVTtoV:= quo<VT|[VT1]>;      // V=T/R is not a ring: multiplication is not well-defined.
                                           // To multiply we need to lift to VT=T/P.
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
            assert2 forall{ v : v in Basis(VT) | mTtoVT((v@@mTtoVT)) eq v };
            assert2 forall{ t : t in ZBasis(T) | t-((mTtoVT(t))@@mTtoVT) in R };

            q:=Index(R,P);
            qpow:=hom<V->V | [mVTtoV(pow(v@@mVTtoV,q)) : v in Basis(V)]>; // This map V->V sends x:->x^q. It is well-def.
            eigen_vals:=[e[1] : e in Setseq(Eigenvalues(Matrix(qpow)))];
            eigen_spaces:=[Kernel(hom<V->V | [qpow(v)-e*v : v in Basis(V)]>): e in eigen_vals]; 

            F:=BaseField(V);
            type1_ws:=[V|]; // the images in V of the generators of the minimal overorders S with dim_F(S/R)=1.
            for E in eigen_spaces do
                n:=Dimension(E);
                // The following loop generates all generators of 1-dimensional subspaces of E.
                for pivot := 1 to n do
                    // Pre-pivot elements are 0, pivot is 1, post-pivot elements vary over F
                    tail_tuples := CartesianPower(F, n - pivot);
                    for t in tail_tuples do
                        w := V!E!([F!0 : j in [1..pivot-1]] cat [F!1] cat [F!x : x in t]);
                        if q eq 2 then
                            // for q eq 2 being a subspace of the eigenspace guarantees that it is mult closed
                            S:=Order([w@@mTV] cat zbR : CheckIsKnownOrder:=CheckIsKnownOrder);
                            Include(~output,S);// necessarily minimal, because it has dim 1
                            Append(~type1_ws,w);
                        else
                            wVT:=w@@mVTtoV;
                            w_square:=pow(wVT,2)@mVTtoV;
                            if Rank(Matrix([w,w_square])) eq 1 then
                                S:=Order(Append(zbR,wVT@@mTtoVT) : CheckIsKnownOrder:=CheckIsKnownOrder);
                                Include(~output,S); // necessarily minimal, because it has dim 1
                                Append(~type1_ws,w);
                            end if;
                        end if;
                    end for;
                end for;
            end for;
            ////////////////////////////////////////////////////////////////////////////
            // The other minimal overorders S of R are such that S/P is a finite field extension of prime degree
            // l>2 of F=R/P. Given a subspace W of V of dimension l-1 whose preimage in VT=T/P is multiplicatively
            // closed, say B=S/P, we need to check that B is a field.
            // If B is not a field, then B contains a non-zero element e which is either nilpotent (with e^2=0) or
            // idempotent (e^2=e). In both cases, F+Fe is a 2-dimensional subalgebra of B, that is, the reduction
            // of a minimal overorder of type 1 found above. Conversely, if B is a field then it contains no
            // proper subalgebra strictly containing F, since [B:F] is prime.
            // Hence S is minimal if and only if W contains none of the lines spanned by type1_ws.
            dims := PrimesUpTo(dV+1); //the plus one is to prevent issues when d=2.
            subs_type2 := Submodules(V : CodimensionLimit := dV-2); //we exclude dim 0 and 1
            for W in subs_type2 do
                if Dimension(W)+1 in dims and
                        not exists{ w : w in type1_ws | w in W } and
                        is_ord_inV(W) then
                    assert Dimension(W) ge 2;
                    S:=Order([v@@mTV:v in Basis(W)] cat zbR : CheckIsKnownOrder:=CheckIsKnownOrder);
                    Include(~output,S);
                end if;
            end for;
        end if;
        assert2 forall{S : S in output | ColonIdeal(R,R!!OneIdeal(S)) eq P};
        // We store the output only if its orders have been checked against the database of known orders.
        // Otherwise, the attribute would contain copies of orders which do not share their attributes
        // with the known ones.
        if CheckIsKnownOrder then
            if not assigned R`MinimalOverOrders then
                R`MinimalOverOrders:=[<P,output>];
            else
                Append(~R`MinimalOverOrders,<P,output>);
            end if;
        end if;
    end if;
    return output;
end intrinsic;


/// Computes the minimal overorders of the order $R$.
intrinsic MinimalOverOrders(R::AlgEtQOrd) -> SetIndx[AlgEtQOrd]
{Computes the minimal overorders of R.}
    output:={@ Universe({@ R @}) | @}; //empty set
    if IsMaximal(R) then
        return output;
    end if;
    // Note: every minimal overorder is a P-MinimalOverOrder for a unique singular prime P.
    // i.e. if P ne Q, then the set of P-minimal overorders and Q-minimal overorders are disjoint!
    return &join[ MinimalOverOrdersAtPrime(R,P):P in SingularPrimes(R) ];
end intrinsic;


// Given an order $R$ and prime $P$ of $R$, returns $R$ and the overorders $S$ of $R$ with conductor $(R:S)$ which is $P$-primary. We recursively produce the minimal $Q$-overorders where $Q$ are primes above $P$. The output is sorted by decreasing `Index(S)`, which is the same as increasing index $[S:R]$, with ties broken using `myHash(S)`. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana.
intrinsic OverOrdersAtPrime(R::AlgEtQOrd, P::AlgEtQIdl) -> SeqEnum[AlgEtQOrd]
{Given an order R and prime P of R, returns R and the overorders S of R with conductor (R:S) which is P-primary. We recursively produce the minimal PP-overorders where PP are primes above P. The output is sorted by decreasing Index(S), which is the same as increasing index [S:R], with ties broken using myHash(S). Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana.}
    require IsPrime(P) : "The ideal P must be prime.";
    require Order(P) eq R : "P must be an ideal of R";

    if assigned R`OverOrdersAtPrimes and exists(output){ tup : tup in R`OverOrdersAtPrimes | tup[1] eq P } then
    // early exit if already computed
        output := output[2];
        return output;
    end if;

    if assigned R`OverOrders then
    // early exit if already computed. This is useful when Loading the data using LoadWKICM
        output:=[ R ] cat [ S : S in R`OverOrders | PrimesAbove(ColonIdeal(R,R!!OneIdeal(S))) eq [ P ] ];
        sort_overorders(~output);
        return output;
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
                pot_new join:=MinimalOverOrdersAtPrime(T,Q);
            end for;
        end for;
        output join:=pot_new;
        done join:=queue;
        queue := pot_new diff done;
    end while;
    output:=Setseq(output);
    sort_overorders(~output);
    assert2 forall{S : S in output | S eq R or PrimesAbove(ColonIdeal(R,R!!OneIdeal(S))) eq [ P ]};
    if not assigned R`OverOrdersAtPrimes then
        R`OverOrdersAtPrimes:=[<P,output>];
    else
        Append(~R`OverOrdersAtPrimes,<P,output>);
    end if;
    return output;
end intrinsic;

/// Returns the overorders of $R$, sorted by decreasing `Index(S)`, which is the same as increasing index $[S:R]$, with ties broken using `myHash`. The parameter `populateoo_in_oo` (default false) determines whether the algorithms fills the attribute `OverOrders` for each overorder of $R$. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana.
intrinsic OverOrders(R::AlgEtQOrd : populateoo_in_oo:=false) -> SeqEnum[AlgEtQOrd]
{Returns the overorders of R, sorted by decreasing Index(S), which is the same as increasing index [S:R], with ties broken using myHash(S). The parameter "populateoo_in_oo" (default false) determines whether we should fill the attribute T`OverOrders for every overorder T of R. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana.}
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
            sort_overorders(~output);
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
{Returns the overorders of R, sorted by decreasing Index(S), which is the same as increasing index [S:R], with ties broken using myHash(S). Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana. The parameter "populateoo_in_oo" (default false) determines whether we should fill the attribute T`OverOrders for every overorder T of R.}
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
    //AttachSpec("../spec");

    SetVerbose("OverOrders",1);
    SetAssertions(2);

    _<x>:=PolynomialRing(Integers());
    f:=(x^4+16);
    A:=EtaleAlgebra(f);
    O:=MaximalOrder(A);
    assert #FindOverOrders(O) eq 1;
    assert #MinimalOverOrders(O) eq 0;

    E:=EquationOrder(A);
    oo:=FindOverOrders(E);
    assert #oo eq 11;

    SetAssertions(1);

    printf " all good!";
*/
