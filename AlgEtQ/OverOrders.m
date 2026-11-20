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
                               OverOrders,        // all overorders
                               TuplesOfOverOrders; // a tuple <oo_at_Ps, tups, pos>, see overorders_as_tuples

sort_key:=func< S | [1/Index(S)] cat myHash(S) >;
// The key used to sort the overorders, see sort_overorders. Distinct orders have distinct keys.

sort_overorders:=procedure(~oo)
// Given a sequence oo of orders, sorts it by decreasing Index(S), breaking ties using myHash(S).
// If the orders in oo are overorders of a fixed order R, then Index(S,R)=Index(R)/Index(S), so this is
// the same as sorting by increasing index [S:R]. In particular, R is always the first entry, and the
// maximal order the last one.
// This makes the output of OverOrdersAtPrime and OverOrders independent of the way the orders are produced.
    keys:=[ sort_key(S) : S in oo ];
    ParallelSort(~keys,~oo);
end procedure;

sum_of_local_overorders:=function(R,oo_at_Ps,c)
// Given an order R, the sequence oo_at_Ps of OverOrdersAtPrime(R,P) for P in SingularPrimes(R), and a sequence c of
// positions, one for each singular prime, returns the order which is the sum of the orders oo_at_Ps[k][c[k]].
// Recall that oo_at_Ps[k][1] eq R. Using this info we can make the generation process a bit faster.
    not_R:=[ k : k in [1..#c] | c[k] ne 1 ]; // the primes where the sum is not equal to R
    if #not_R eq 0 then
        return R;
    elif #not_R eq 1 then
        // this case happens for the minimal overorders of R, and when there is only one singular prime of R
        k:=not_R[1];
        return oo_at_Ps[k][c[k]];
    else
        // The orders are overorders of R at distinct primes of R: their sum is an order, which is the Z-span of
        // their Z-bases. Hence we can use Check:=0.
        return Order(&cat[ ZBasis(oo_at_Ps[k][c[k]]) : k in not_R ] : Check:=0);
    end if;
end function;

overorders_as_tuples:=function(R)
// Given a non-maximal order R with singular primes P_1,...,P_n, returns:
//  - oo_at_Ps: the sequence of OverOrdersAtPrime(R,P_k) for k in [1..n];
//  - tups: a sequence of sequences of integers, such that OverOrders(R)[i] is the sum of the orders
//          oo_at_Ps[k][tups[i][k]], for k in [1..n];
//  - pos: an AssociativeArray such that pos[tups[i]] eq i.
// This describes the bijection L_1 x ... x L_n -> OverOrders(R), where L_k is the set of P_k-overorders of R and R
// itself, given by the sum of the orders: see OverOrders. It respects the inclusions in both directions.
// The output is stored in the attribute R`TuplesOfOverOrders. It is filled by OverOrders when it computes the
// overorders. If the overorders have been assigned in another way, for example by WKICM or LoadWKICM, we rebuild the
// order corresponding to each tuple (one HNF each) and we find its position in OverOrders(R) with Index.
    if not assigned R`TuplesOfOverOrders then
        oo:=OverOrders(R); // it might assign R`TuplesOfOverOrders
        if not assigned R`TuplesOfOverOrders then
            pp:=SingularPrimes(R);
            n:=#pp;
            oo_at_Ps:=[ OverOrdersAtPrime(R,P) : P in pp ];
            assert forall{ x : x in oo_at_Ps | x[1] eq R };
            oo_set:={@ S : S in oo @}; // an indexed set, to find the position of an order using its Hash
            tups:=[ [] : i in [1..#oo] ];
            pos:=AssociativeArray();
            for c0 in CartesianProduct([ [1..#x] : x in oo_at_Ps ]) do
                c:=[ c0[k] : k in [1..n] ];
                i:=Index(oo_set,sum_of_local_overorders(R,oo_at_Ps,c));
                assert i ne 0;
                tups[i]:=c;
                pos[c]:=i;
            end for;
            assert #Keys(pos) eq #oo; // every overorder corresponds to a tuple
            R`TuplesOfOverOrders:=<oo_at_Ps,tups,pos>;
        end if;
    end if;
    return Explode(R`TuplesOfOverOrders);
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
    // Early exit if all the overorders are already known, for example when they have been assigned by WKICM or
    // LoadWKICM. We get the P-overorders by selecting them from R`OverOrders.
    // An overorder S of R is a P-overorder if and only if the R-module S/R is supported only at P. The support of S/R
    // consists of the primes containing its annihilator, which is the colon ideal (R:S). These primes are singular,
    // since S/R is zero locally at the other primes.
    // We sort all the overorders according to the support of S/R in one pass, and we store the result for every
    // singular prime of R: each overorder is tested only once, whichever prime is asked next.
        pp:=SingularPrimes(R);
        iP:=Index(pp,P);
        if iP eq 0 then
            return [ R ]; // R is maximal at P
        end if;
        oo_at_Ps:=[ [ R ] : Q in pp ];
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
        for k->Q in pp do
            output:=oo_at_Ps[k];
            sort_overorders(~output);
            oo_at_Ps[k]:=output;
            if not assigned R`OverOrdersAtPrimes then
                R`OverOrdersAtPrimes:=[<Q,output>];
            elif not exists{ tup : tup in R`OverOrdersAtPrimes | tup[1] eq Q } then
                Append(~R`OverOrdersAtPrimes,<Q,output>);
            end if;
        end for;
        return oo_at_Ps[iP];
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
            // Let P_1,...,P_n be the singular primes of R. For each k, let L_k be the set consisting of R and of
            // the P_k-overorders of R, that is, the overorders S of R such that S/R is supported only at P_k.
            // The map
            //      L_1 x ... x L_n -> { overorders of R }       (S_1,...,S_n) :-> S_1+...+S_n
            // is a bijection which respects inclusions in both directions: S_1+...+S_n is the unique overorder
            // of R which is equal to S_k locally at P_k, for each k.
            // We encode an overorder by the tuple of the positions of the S_k's in the sequences oo_at_Ps[k],
            // and we store the correspondence in the attribute R`TuplesOfOverOrders (see overorders_as_tuples).
            pp:=SingularPrimes(R);
            n:=#pp;
            vtime OverOrders,2 : oo_at_Ps:=[ OverOrdersAtPrime(R,P) : P in pp ];
            vprintf OverOrders,2 : "P-overorders for each P %o\n",[#x : x in oo_at_Ps];
            vprintf OverOrders,2 : "total number of overorders %o\n",&*[#x : x in oo_at_Ps];
            assert forall{x : x in oo_at_Ps | x[1] eq R}; // used by sum_of_local_overorders
            tups:=[ [ c[k] : k in [1..n] ] : c in CartesianProduct([[1..#x] : x in oo_at_Ps]) ];
            output:=[ sum_of_local_overorders(R,oo_at_Ps,c) : c in tups ];
            // we sort output as in sort_overorders, and tups in the same way
            keys:=[ sort_key(S) : S in output ];
            perm:=[1..#output];
            ParallelSort(~keys,~perm);
            output:=[ output[i] : i in perm ];
            tups:=[ tups[i] : i in perm ];
            pos:=AssociativeArray();
            for i->c in tups do
                pos[c]:=i;
            end for;
            assert2 R in output;
            assert2 MaximalOrder(Algebra(R)) in output;
            R`TuplesOfOverOrders:=<oo_at_Ps,tups,pos>;
        end if;
        R`OverOrders:=output;
    end if;
    if populateoo_in_oo and not IsMaximal(R) then
        // For S and T in L_1 x ... x L_n, we have S subset T if and only if S_k subset T_k for every k, since the
        // inclusion can be checked locally at each P_k, where S and T are equal to S_k and T_k.
        // Hence we only need the inclusions in each L_k: #L_1^2+...+#L_n^2 inclusion tests instead of
        // (#L_1*...*#L_n)^2.
        oo:=R`OverOrders;
        oo_at_Ps,tups,pos:=overorders_as_tuples(R);
        n:=#oo_at_Ps;
        // up[k][a] is the sequence of the positions b such that oo_at_Ps[k][a] subset oo_at_Ps[k][b], a included.
        up:=[ [ [ b : b in [1..#L] | L[a] subset L[b] ] : a in [1..#L] ] : L in oo_at_Ps ];
        for i in [1..#oo] do
            S:=oo[i];
            if not assigned S`OverOrders then
                c:=tups[i];
                // the positions of the overorders of S, sorted as in OverOrders(R), hence as in OverOrders(S)
                inds:=Sort([ pos[[ d[k] : k in [1..n] ]] : d in CartesianProduct([ up[k][c[k]] : k in [1..n] ]) ]);
                S`OverOrders:=[ oo[j] : j in inds ];
                assert3 S`OverOrders eq [ T : T in oo | S subset T ]; // quadratic in #oo
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

    // populateoo_in_oo, computed from the overorders at each singular prime, against the inclusion test between all
    // pairs of overorders; and OverOrdersAtPrime when the overorders are already known (after LoadWKICM), against the
    // computation from scratch. With 1, 2 and 4 singular primes.
    for f in [ x^4+16, x^4+4*x^3+6*x^2+44*x+121, (x^2+5)*(x^2+20) ] do
        A:=EtaleAlgebra(f);
        E:=EquationOrder(A);
        oo:=OverOrders(E : populateoo_in_oo:=true);
        assert forall{ S : S in oo | S`OverOrders eq [ T : T in oo | S subset T ] };
        pp:=SingularPrimes(E);
        oo_at_Ps:=[ OverOrdersAtPrime(E,P) : P in pp ];
        E2:=LoadWKICM(PrintWKICM(E)); // a new algebra, where only the overorders and the weak classes are known
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
    end for;

    SetAssertions(1);

    printf " all good!";
*/
