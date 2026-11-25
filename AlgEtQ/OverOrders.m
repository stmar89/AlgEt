/////////////////////////////////////////////////////
// Copyright 2025.
// Stefano Marseglia, stefano.marseglia89@gmail.com
// https://stmar89.github.io/index.html
//
// Distributed under the terms of the CC-BY 4.0 licence.
// https://creativecommons.org/licenses/by/4.0/
/////////////////////////////////////////////////////


freeze;

declare verbose OverOrders, 1;

import "Ord.m" : op_over_product;

// Attributes of an order R, with the only functions that write them.
declare attributes AlgEtQOrd : MinimalOverOrders, // An AssociativeArray keyed by the singular
                                                  // primes P of R: the value at P is the indexed
                                                  // set of all the minimal overorders T of R with
                                                  // conductor (R:T) = P. Keyed by the ideal, using
                                                  // the Hash of ideals: P need not be the same
                                                  // object as in SingularPrimes(R). Written by
                                                  // MinimalOverOrdersAtPrime.
                               LocalOverOrders,   // A record in local_overorders_fmt: the
                                                  // overorders of R at each singular prime, and
                                                  // how their sums correspond to OverOrders(R).
                                                  // Field at_primes written by OverOrdersAtPrime,
                                                  // field pos by OverOrders and
                                                  // positions_of_overorders.
                               OverOrders;        // All the overorders of R, sorted with
                                                  // sort_overorders. Written by OverOrders (also
                                                  // for the overorders of R, with
                                                  // populateoo_in_oo) and by store_overorders,
                                                  // which WKICM and LoadWKICM use to store the
                                                  // overorders they find.

local_overorders_fmt:=recformat< at_primes, pos >;
// Let P_1,...,P_n be the singular primes of R, in the order of SingularPrimes(R), and let L_k be
// the set consisting of R and of the P_k-overorders of R, that is, the overorders S of R such that
// S/R is supported only at P_k.
//  - at_primes: an AssociativeArray: at_primes[k] is L_k as a sequence, sorted with
//               sort_overorders (so R is its first entry). The entry k is defined once
//               OverOrdersAtPrime(R,P_k) has been computed: the primes are computed one at a time,
//               unless the overorders of R are already known.
//  - pos:       an AssociativeArray, assigned once all the L_k's and OverOrders(R) are known: for
//               every sequence c=[c_1,...,c_n] with c_k in [1..#L_k], OverOrders(R)[pos[c]] is the
//               sum of the orders L_k[c_k].
// By Hofmann-Sircana, the map
//      L_1 x ... x L_n -> { overorders of R }       (S_1,...,S_n) :-> S_1+...+S_n
// is a bijection which respects inclusions in both directions: S_1+...+S_n is the unique overorder
// of R which is equal to S_k locally at P_k, for each k. The overorders are encoded by positions in
// the L_k's, rather than by the sequences of orders [S_1,...,S_n], because looking up a sequence of
// integers in an AssociativeArray is about 5 times faster than looking up a sequence of orders, and
// populateoo_in_oo does millions of lookups.

sort_key:=func< S | [1/Index(S)] cat myHash(S) >;
// The key used to sort the overorders, see sort_overorders. Distinct orders have distinct keys.

sort_overorders:=procedure(~oo)
// Given a sequence oo of orders, sorts it by decreasing Index(S), breaking ties using myHash(S).
// If the orders in oo are overorders of a fixed order R, then Index(S,R)=Index(R)/Index(S), so this
// is the same as sorting by increasing index [S:R]. In particular, R is always the first entry, and
// the maximal order the last one.
// This makes the output of OverOrdersAtPrime and OverOrders independent of the way the orders are
// produced.
    keys:=[ sort_key(S) : S in oo ];
    ParallelSort(~keys,~oo);
end procedure;

store_overorders:=procedure(R,oo)
// Given an order R and the sequence oo of all its overorders, in any order, stores them in
// R`OverOrders, sorted with sort_overorders. If R`OverOrders is already assigned, it is kept: the
// objects in it may carry attributes, and the positions in R`LocalOverOrders`pos refer to its
// order. This is how the overorders found by other means (WKICM, LoadWKICM) are stored.
    if not assigned R`OverOrders then
        sort_overorders(~oo);
        assert2 R eq oo[1];
        R`OverOrders:=oo;
    else
        assert2 #oo eq #R`OverOrders and { S : S in oo } eq { S : S in R`OverOrders };
    end if;
end procedure;

local_overorders:=function(R)
// Returns the record R`LocalOverOrders, after creating it empty if needed.
    if not assigned R`LocalOverOrders then
        R`LocalOverOrders:=rec< local_overorders_fmt | at_primes:=AssociativeArray() >;
    end if;
    return R`LocalOverOrders;
end function;

p_overorders_from_known:=function(R)
// Given a non-maximal order R whose overorders are known (R`OverOrders is assigned), returns the
// sequence of the L_k's, see local_overorders_fmt, selected from R`OverOrders.
// An overorder S of R is a P-overorder if and only if the R-module S/R is supported only at P. The
// support of S/R consists of the primes containing its annihilator, which is the colon ideal (R:S).
// These primes are singular, since S/R is zero locally at the other primes. So each overorder is
// tested once, for all the singular primes.
    pp:=SingularPrimes(R);
    oo_at_Ps:=[ [ R ] : P in pp ];
    for S in R`OverOrders do
        if S ne R then
            cS:=ColonIdeal(R,R!!OneIdeal(S));
            supp:=[ k : k in [1..#pp] | cS subset pp[k] ];
            assert #supp ge 1; // since S ne R
            if #supp eq 1 then
                Append(~oo_at_Ps[supp[1]],S);
            end if;
        end if;
    end for;
    for k in [1..#pp] do
        L:=oo_at_Ps[k];
        sort_overorders(~L);
        oo_at_Ps[k]:=L;
    end for;
    return oo_at_Ps;
end function;

p_overorders_from_minimal:=function(R,P)
// Given an order R and a singular prime P of R, returns R and the P-overorders of R, sorted with
// sort_overorders. They are produced recursively from R, adding at each step the minimal
// Q-overorders of the orders found so far, for the primes Q above P (Hofmann-Sircana).
    ppO:=PrimesAbove(MaximalOrder(Algebra(R))!!P);
    queue := {@ R @};
    output:={@ R @};
    done:={@ @};
    while #queue gt 0 do
        pot_new:={@ @};
        for T in queue do
            pp:={@ PrimeBelow(T,Q) : Q in ppO @}; // the primes of T above P
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
    return output;
end function;

sum_of_two_overorders:=function(R,S,T)
// Given overorders S and T of R, which differ from R at disjoint sets of singular primes, returns
// S+T: an order, the Z-span of their Z-bases, since at each prime it is locally equal to S or to T.
// If one of them is R (the first entry of each L_k, see local_overorders_fmt), the sum is the other
// one.
    if T eq R then
        return S;
    elif S eq R then
        return T;
    end if;
    return Order(ZBasis(S) cat ZBasis(T) : Check:=0);
end function;

all_sums_of_local_overorders:=function(R,oo_at_Ps)
// Returns the sequence cs of all the sequences of positions c=[c_1,...,c_n] with c_k in
// [1..#oo_at_Ps[k]], and the sequence sums, where sums[i] is the sum of the orders
// oo_at_Ps[k][cs[i][k]]. The sums are computed with shared prefixes (see op_over_product): the sum
// for c is the sum for the prefix c_1,...,c_(n-1), which is computed once, plus oo_at_Ps[n][c_n];
// each sum is one HNF of 2*dim elements, instead of n*dim. When the entries c_k are 1, the sum is
// R, or a P-overorder, which is returned as it is.
    n:=#oo_at_Ps;
    cs:=[ [ c[k] : k in [1..n] ] : c in CartesianProduct([ [1..#x] : x in oo_at_Ps ]) ];
    sums:=op_over_product(oo_at_Ps,func< S,T | sum_of_two_overorders(R,S,T) >);
    assert #sums eq #cs;
    return cs,sums;
end function;

positions_of_overorders:=function(R)
// Given a non-maximal order R, returns the sequence oo_at_Ps of the L_k's and the AssociativeArray
// pos, see local_overorders_fmt. If OverOrders(R) computed the overorders, pos is already known.
// Otherwise, for example if the overorders were stored by WKICM or LoadWKICM, every sum is computed
// again (one HNF each) and looked up in OverOrders(R) with an indexed set, which uses the Hash of
// the orders. The result is stored in R`LocalOverOrders.
    oo:=OverOrders(R);
    oo_at_Ps:=[ OverOrdersAtPrime(R,P) : P in SingularPrimes(R) ];
    data:=local_overorders(R);
    if not assigned data`pos then
        cs,sums:=all_sums_of_local_overorders(R,oo_at_Ps);
        oo_set:={@ S : S in oo @};
        pos:=AssociativeArray();
        for i->c in cs do
            j:=Index(oo_set,sums[i]);
            assert j ne 0;
            pos[c]:=j;
        end for;
        assert #oo_set eq #cs; // the bijection
        data`pos:=pos;
        R`LocalOverOrders:=data;
    end if;
    return oo_at_Ps,data`pos;
end function;

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

    if assigned R`MinimalOverOrders and IsDefined(R`MinimalOverOrders,P) then
        return R`MinimalOverOrders[P]; // already computed
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
            V,mVTtoV:= quo<VT|[VT1]>; // V=T/R is not a ring: multiplication is not well-defined.
                                      // To multiply we need to lift to VT=T/P.
            dV := Dimension(V);
            assert #dVT eq dV+1; 
            // minimal P-overorder S of R such that S/R has dimension 1 over R/P satisfy the
            // following property: S/R is necessarily contained in an eigenspace of x->x^q acting on
            // V=T/R.
            mTV := map<T->V|x:->mVTtoV(mTtoVT(x)),y:->y@@mVTtoV@@mTtoVT>; // mTV: T->T/R=V
            // We create 3 auxiliary function: a multiplication map on VT, a n-power map on VT, and
            // a functions that checks whether a subspace of V is mult closed, that is, the
            // reduction of an order. The latter builds upon the first 2.
            // This approach keeps the coefficients bounded, allowing to use linear algebra over
            // finite fields rather than over the integers.
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
                // Given a subvectors space W of V, returns whether its preimage in VT is
                // multiplicatively closed. This happens precisely when W is the reduction of an
                // order S such that R < S < T.
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
            // This map V->V sends x:->x^q. It is well-def.
            qpow:=hom<V->V | [mVTtoV(pow(v@@mVTtoV,q)) : v in Basis(V)]>;
            eigen_vals:=[e[1] : e in Setseq(Eigenvalues(Matrix(qpow)))];
            eigen_spaces:=[Kernel(hom<V->V | [qpow(v)-e*v : v in Basis(V)]>): e in eigen_vals]; 

            F:=BaseField(V);
            type1_ws:=[V|]; // the images in V of the generators of the minimal overorders S with
                            // dim_F(S/R)=1.
            for E in eigen_spaces do
                n:=Dimension(E);
                // The following loop generates all generators of 1-dimensional subspaces of E.
                for pivot := 1 to n do
                    // Pre-pivot elements are 0, pivot is 1, post-pivot elements vary over F
                    tail_tuples := CartesianPower(F, n - pivot);
                    for t in tail_tuples do
                        w := V!E!([F!0 : j in [1..pivot-1]] cat [F!1] cat [F!x : x in t]);
                        if q eq 2 then
                            // for q eq 2 being a subspace of the eigenspace guarantees that it
                            // is mult closed
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
            // The other minimal overorders S of R are such that S/P is a finite field extension of
            // prime degree l>2 of F=R/P. Given a subspace W of V of dimension l-1 whose preimage in
            // VT=T/P is multiplicatively closed, say B=S/P, we need to check that B is a field.
            // If B is not a field, then B contains a non-zero element e which is either nilpotent
            // (with e^2=0) or idempotent (e^2=e). In both cases, F+Fe is a 2-dimensional subalgebra
            // of B, that is, the reduction of a minimal overorder of type 1 found above.
            // Conversely, if B is a field then it contains no proper subalgebra strictly containing
            // F, since [B:F] is prime.
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
        // We store the output only if its orders have been checked against the database of known
        // orders. Otherwise, the attribute would contain copies of orders which do not share their
        // attributes with the known ones.
        if CheckIsKnownOrder then
            mins:=assigned R`MinimalOverOrders select R`MinimalOverOrders else AssociativeArray();
            mins[P]:=output;
            R`MinimalOverOrders:=mins;
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


/// Given an order $R$ and prime $P$ of $R$, returns $R$ and the overorders $S$ of $R$ with conductor $(R:S)$ which is $P$-primary. We recursively produce the minimal $Q$-overorders where $Q$ are primes above $P$. The output is sorted by decreasing `Index(S)`, which is the same as increasing index $[S:R]$, with ties broken using `myHash(S)`. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana.
intrinsic OverOrdersAtPrime(R::AlgEtQOrd, P::AlgEtQIdl) -> SeqEnum[AlgEtQOrd]
{Given an order R and prime P of R, returns R and the overorders S of R with conductor (R:S) which is P-primary. We recursively produce the minimal PP-overorders where PP are primes above P. The output is sorted by decreasing Index(S), which is the same as increasing index [S:R], with ties broken using myHash(S). Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana.}
    require IsPrime(P) : "The ideal P must be prime.";
    require Order(P) eq R : "P must be an ideal of R";
    k:=Index(SingularPrimes(R),P);
    if k eq 0 then
        return [ R ]; // R is maximal at P
    end if;
    data:=local_overorders(R);
    if not IsDefined(data`at_primes,k) then
        if assigned R`OverOrders then
            // The overorders are known, for example from WKICM or LoadWKICM: we select the
            // P-overorders from them, for all the singular primes at once.
            for j->L in p_overorders_from_known(R) do
                if not IsDefined(data`at_primes,j) then
                    data`at_primes[j]:=L;
                end if;
            end for;
        else
            data`at_primes[k]:=p_overorders_from_minimal(R,P);
        end if;
        R`LocalOverOrders:=data;
    end if;
    return data`at_primes[k];
end intrinsic;

/// Returns the overorders of $R$, sorted by decreasing `Index(S)`, which is the same as increasing index $[S:R]$, with ties broken using `myHash`. The parameter `populateoo_in_oo` (default false) determines whether the algorithms fills the attribute `OverOrders` for each overorder of $R$. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana.
intrinsic OverOrders(R::AlgEtQOrd : populateoo_in_oo:=false) -> SeqEnum[AlgEtQOrd]
{Returns the overorders of R, sorted by decreasing Index(S), which is the same as increasing index [S:R], with ties broken using myHash(S). The parameter "populateoo_in_oo" (default false) determines whether we should fill the attribute T`OverOrders for every overorder T of R. Based on "On the computations of overorders" by Tommy Hofmann and Carlo Sircana.}
    if not assigned R`OverOrders then
        if IsMaximal(R) then
            R`OverOrders:=[R];
        else
            // The overorders of R are the sums of its overorders at each singular prime: see
            // local_overorders_fmt.
            oo_at_Ps:=[ OverOrdersAtPrime(R,P) : P in SingularPrimes(R) ];
            vprintf OverOrders,2 : "P-overorders for each P %o\n",[#x : x in oo_at_Ps];
            vprintf OverOrders,2 : "total number of overorders %o\n",&*[#x : x in oo_at_Ps];
            // all_sums_of_local_overorders uses that R is the first entry of each L_k
            assert forall{ x : x in oo_at_Ps | x[1] eq R };
            cs,sums:=all_sums_of_local_overorders(R,oo_at_Ps);
            // We sort the sums as in sort_overorders, and cs in the same way: then their positions
            // are known.
            keys:=[ sort_key(S) : S in sums ];
            perm:=[1..#sums];
            ParallelSort(~keys,~perm); // the i-th overorder is sums[perm[i]]
            pos:=AssociativeArray();
            for i->j in perm do
                pos[cs[j]]:=i;
            end for;
            oo:=[ sums[j] : j in perm ];
            assert2 R in oo;
            assert2 MaximalOrder(Algebra(R)) in oo;
            R`OverOrders:=oo;
            data:=local_overorders(R);
            data`pos:=pos;
            R`LocalOverOrders:=data;
        end if;
    end if;
    if populateoo_in_oo and not IsMaximal(R) then
        // Write S=S_1+...+S_n and T=T_1+...+T_n as in local_overorders_fmt. Then S subset T if and
        // only if S_k subset T_k for every k, since the inclusion can be checked locally at each
        // P_k, where S and T are equal to S_k and T_k. Hence we only need the inclusions in each
        // L_k: #L_1^2+...+#L_n^2 inclusion tests instead of (#L_1*...*#L_n)^2.
        oo:=R`OverOrders;
        oo_at_Ps,pos:=positions_of_overorders(R);
        n:=#oo_at_Ps;
        // up[k][a] is the sequence of the positions b such that oo_at_Ps[k][a] subset
        // oo_at_Ps[k][b], a included.
        up:=[ [ [ b : b in [1..#L] | L[a] subset L[b] ] : a in [1..#L] ] : L in oo_at_Ps ];
        for c in Keys(pos) do
            S:=oo[pos[c]];
            if not assigned S`OverOrders then
                // the positions of the overorders of S, sorted as in OverOrders(R), hence as in
                // OverOrders(S)
                inds:=Sort([ pos[[ d[k] : k in [1..n] ]] : d in CartesianProduct([ up[k][c[k]] : k in [1..n] ]) ]);
                S`OverOrders:=[ oo[j] : j in inds ];
                assert3 S`OverOrders eq [ T : T in oo | S subset T ]; // quadratic in #oo
            end if;
        end for;
    end if;
    return R`OverOrders;
end intrinsic;

///ditto
intrinsic FindOverOrders(R::AlgEtQOrd : populateoo_in_oo:=false) -> SeqEnum[AlgEtQOrd]
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
    // the minimal overorders are stored keyed by the prime: an equal prime which is a different
    // object finds them
    P:=SingularPrimes(E)[1];
    min:=MinimalOverOrdersAtPrime(E,P);
    P2:=Ideal(E,ZBasis(P));
    assert not IsIdentical(P,P2) and P2 eq P;
    assert IsDefined(E`MinimalOverOrders,P2) and MinimalOverOrdersAtPrime(E,P2) eq min;

    // populateoo_in_oo, computed from the overorders at each singular prime, against the inclusion
    // test between all pairs of overorders; and OverOrdersAtPrime when the overorders are already
    // known (after LoadWKICM), against the computation from scratch. With 1, 2 and 4 singular
    // primes.
    for f in [ x^4+16, x^4+4*x^3+6*x^2+44*x+121, (x^2+5)*(x^2+20) ] do
        A:=EtaleAlgebra(f);
        E:=EquationOrder(A);
        oo:=OverOrders(E : populateoo_in_oo:=true);
        assert forall{ S : S in oo | S`OverOrders eq [ T : T in oo | S subset T ] };
        pp:=SingularPrimes(E);
        oo_at_Ps:=[ OverOrdersAtPrime(E,P) : P in pp ];
        E2:=LoadWKICM(PrintWKICM(E)); // a new algebra, where only the overorders and the weak
                                      // classes are known
        pp2:=SingularPrimes(E2);
        assert #pp2 eq #pp;
        for P2 in pp2 do
            // the primes of E2 match the ones of E through the Z-coordinates of their bases
            i:=[ i : i in [1..#pp] | myHash(pp[i]) eq myHash(P2) ];
            assert #i eq 1;
            assert [ myHash(S) : S in OverOrdersAtPrime(E2,P2) ] eq [ myHash(S) : S in oo_at_Ps[i[1]] ];
        end for;
        oo2:=OverOrders(E2 : populateoo_in_oo:=true);
        assert forall{ S : S in oo2 | S`OverOrders eq [ T : T in oo2 | S subset T ] };
        // OverOrdersAtPrime at one prime only, then WKICM, which stores the overorders: the
        // overorders at the other primes are then selected from them. Everything must be as
        // computed directly.
        E3:=EquationOrder(EtaleAlgebra(f));
        pp3:=SingularPrimes(E3);
        _:=OverOrdersAtPrime(E3,pp3[1]);
        _:=WKICM(E3);
        for P3 in pp3 do
            i:=[ i : i in [1..#pp] | myHash(pp[i]) eq myHash(P3) ];
            assert #i eq 1;
            assert [ myHash(S) : S in OverOrdersAtPrime(E3,P3) ] eq [ myHash(S) : S in oo_at_Ps[i[1]] ];
        end for;
        assert [ myHash(S) : S in OverOrders(E3) ] eq [ myHash(S) : S in oo ];
    end for;

    SetAssertions(1);

    printf " all good!";
*/
