/////////////////////////////////////////////////////
// Copyright 2025.
// Stefano Marseglia, stefano.marseglia89@gmail.com
// https://stmar89.github.io/index.html
//
// Distributed under the terms of the CC-BY 4.0 licence.
// https://creativecommons.org/licenses/by/4.0/
/////////////////////////////////////////////////////



freeze;

declare verbose GraphOverOrders, 1;

declare attributes AlgEtQOrd : GraphOverOrders;

import "OverOrders.m" : positions_of_overorders;

///## Lattice of inclusions of the overorders

//------------
// GraphOverOrders : the graph of (minimal) inclusions of overorders of R
//------------

/// Given an order $R$, returns the graph $G$ of minimal inclusions of the overorders of $R$. The vertices of $G$ are integers between $1$ and the number of overorders of $R$, and there is an edge $[i,j]$ if and only if the $j$-th overorder of $R$ is a minimal overorder of the $i$-th.
intrinsic GraphOverOrders(R::AlgEtQOrd) -> GrphDir
{Given an order R, returns the graph G of minimal inclusions of the overorders of R. The vertices of G are integers between 1 and the number of OverOrders(R), and there is an edge [i,j] if and only if OverOrders(R)[j] is a minimal overorder of OverOrders(R)[i].}
    if not assigned R`GraphOverOrders then
        if IsMaximal(R) then
            G:=Digraph< 1 | : SparseRep:=true >; // sparse, as below
        else
            // Let P_1,...,P_n be the singular primes of R. For each k, let L_k be the set
            // consisting of R and of the P_k-overorders of R, that is, the overorders S of R such
            // that S/R is supported only at P_k. By Hofmann-Sircana "On the computation of
            // overorders", the map
            //      L_1 x ... x L_n -> { overorders of R }       (S_1,...,S_n) :-> S_1+...+S_n
            // is a bijection which respects inclusions in both directions: S_1+...+S_n is the
            // unique overorder of R which is equal to S_k locally at P_k, for each k. It follows
            // that T_1+...+T_n is a minimal overorder of S_1+...+S_n if and only if
            //  - S_k eq T_k for all indices k but one, say k0, and
            //  - T_k0 is a minimal overorder of S_k0.
            // Hence we get the graph of R from the graphs of the L_k's, which are much smaller. We
            // proceed in 3 steps:
            //  1) we get the L_k's, and the position in OverOrders(R) of the overorder corresponding
            //     to each tuple;
            //  2) we compute the minimal inclusions in each L_k;
            //  3) we combine the minimal inclusions of the L_k's into the minimal inclusions of the
            //     overorders of R.
            O:=MaximalOrder(Algebra(R));
            pp:=SingularPrimes(R);
            n:=#pp;

            // Step 1: oo_at_Ps[k] is L_k, as a sequence, with R as first entry. An overorder is
            // encoded by a sequence of positions c=[c_1,...,c_n], which stands for the sum of the
            // orders oo_at_Ps[k][c_k], and pos[c] is its position in OverOrders(R). See
            // local_overorders_fmt and positions_of_overorders in OverOrders.m: pos is computed by
            // OverOrders(R), or when it is needed for the first time if the overorders were found
            // in another way, for example by WKICM or LoadWKICM.
            vprintf GraphOverOrders,1 : "Computing the overorders at each singular prime, as sums...";
            t0:=Cputime();
            output:=OverOrders(R);
            oo_at_Ps,pos:=positions_of_overorders(R);
            assert #Keys(pos) eq #output;
            vprintf GraphOverOrders,1 : "Done in %o secs.\n",Cputime(t0);
            vprintf GraphOverOrders,1 : "Number of P-overorders for each singular prime P: %o\n",[ #x : x in oo_at_Ps ];

            // Step 2: min_at_Ps[k][i] is the sequence of the positions j such that oo_at_Ps[k][j]
            // is a minimal overorder of oo_at_Ps[k][i].
            // Let S be in L_k, with P=P_k. If T is a minimal overorder of S then the conductor
            // (S:T) is a prime Q of S, and T is in L_k if and only if Q contains P. So the minimal
            // overorders of S in L_k are the output of MinimalOverOrdersAtPrime(S,Q), where Q runs
            // over the primes of S above P.
            // We do not read them from the attribute S`MinimalOverOrders, because it cannot be
            // trusted to be complete: it is filled as a side effect of OverOrdersAtPrime(R,P) only
            // when this intrinsic computes the P-overorders from scratch, and this does not happen
            // if R`OverOrders is already assigned, for example by WKICM or by LoadWKICM.
            // MinimalOverOrdersAtPrime returns the stored minimal overorders if they are known, and
            // computes (and stores) them otherwise.
            vprintf GraphOverOrders,1 : "Computing the minimal inclusions at each singular prime...";
            t0:=Cputime();
            min_at_Ps:=[];
            for k->P in pp do
                oo_P:=oo_at_Ps[k];
                oo_P_set:={@ S : S in oo_P @}; // an indexed set, to find the position of an order
                                               // using its Hash
                ppO:=PrimesAbove(O!!P); // the primes of O above P
                min_P:=[];
                for S in oo_P do
                    // The primes of S above P are the primes of S below the primes of O above P.
                    // Distinct primes of O might be above the same prime of S: the indexed set
                    // removes the repetitions. PrimeBelow returns the primes stored in S by
                    // OverOrdersAtPrime, so the minimal overorders stored in S under them are found
                    // at once.
                    qq:={@ PrimeBelow(S,QO) : QO in ppO @};
                    min_S:=[];
                    for Q in qq do
                        // The next line returns the empty set if S is maximal at Q.
                        for T in MinimalOverOrdersAtPrime(S,Q) do
                            j:=Index(oo_P_set,T);
                            assert j ne 0; // T is a P-overorder of R
                            Append(~min_S,j);
                        end for;
                    end for;
                    Append(~min_P,min_S);
                end for;
                // tests
                assert2 forall{ i : i in [1..#oo_P] | forall{ j : j in min_P[i] |
                                oo_P[i] subset oo_P[j] } }; // inclusions
                // the next test is very time consuming: cubic in #oo_P
                assert3 forall{ i : i in [1..#oo_P] | forall{ j : j in min_P[i] |
                                [ S : S in oo_P | S ne oo_P[i] and oo_P[i] subset S and S subset oo_P[j] ]
                                    eq [ oo_P[j] ] } }; // minimal inclusions
                Append(~min_at_Ps,min_P);
            end for;
            vprintf GraphOverOrders,1 : "Done in %o secs.\n",Cputime(t0);
            vprintf GraphOverOrders,1 : "Number of minimal inclusions for each singular prime: %o\n",
                                        [ &+[ #x : x in min_P ] : min_P in min_at_Ps ];

            // Step 3: the minimal overorders of the overorder corresponding to the tuple c
            // correspond to the tuples obtained from c by replacing exactly one entry c[k] with an
            // element of min_at_Ps[k][c[k]].
            vprintf GraphOverOrders,1 : "Computing the minimal inclusions of the overorders...";
            t0:=Cputime();
            min_inclusions:=[];
            for c in Keys(pos) do
                i:=pos[c];
                for k in [1..n] do
                    for j in min_at_Ps[k][c[k]] do
                        d:=c;
                        d[k]:=j;
                        Append(~min_inclusions,[i,pos[d]]);
                    end for;
                end for;
            end for;
            vprintf GraphOverOrders,1 : "Done in %o secs.\n",Cputime(t0);
            // tests
            // inclusions
            assert2 forall{ incl : incl in min_inclusions | output[incl[1]] subset output[incl[2]] };
            // minimal inclusions; the next test is very time consuming: cubic in #output
            assert3 forall{ incl : incl in min_inclusions | [ S : S in output |
                            S ne output[incl[1]] and output[incl[1]] subset S
                                and S subset output[incl[2]] ] eq [ output[incl[2]] ] };
            // The sparse representation: the dense one (an adjacency matrix) is refused by Magma
            // for 100000 vertices or more ("Dimension too big to create graph"), and these graphs
            // have only a few edges per vertex.
            G:=Digraph< #output | min_inclusions : SparseRep:=true >;
        end if;
        R`GraphOverOrders:=G;
    end if;
    return R`GraphOverOrders;
end intrinsic;

/* TESTS

    printf "### Testing GraphOverOrders:";
    //AttachSpec("../spec");
    SetAssertions(2);
    _<x>:=PolynomialRing(Integers());

    // the edges of a graph as a set of pairs of integers
    edges:=func< G | { [Index(InitialVertex(e)),Index(TerminalVertex(e))] : e in Edges(G) } >;

    min_inclusions_naive:=function(oo)
    // Given a sequence oo of orders, returns the set of pairs [i,j] such that oo[i] is strictly
    // contained in oo[j] and no other order of oo lies in between. It tests all inclusions:
    // quadratic in #oo.
        above:=[ { j : j in [1..#oo] | j ne i and oo[i] subset oo[j] } : i in [1..#oo] ];
        output:={};
        for i in [1..#oo] do
            for j in above[i] do
                if not exists{ k : k in above[i] | j in above[k] } then
                    Include(~output,[i,j]);
                end if;
            end for;
        end for;
        return output;
    end function;

    fs:=[
          x^8 + 16, //1 sing prime
          x^4-10000*x^3-10000*x^2-10000*x-10000, //2 sing primes
          x^4-30^3*x^3-30^3*x^2-30^3*x-30^3, //3 sing primes
        ];
    data := [<0, 1, 89, 186>, <0, 1, 297, 950>, <0, 1, 64, 144>];
    for i->f in fs do
        printf ".";
        A:=EtaleAlgebra(f);
        tM := GraphOverOrders(MaximalOrder(A));
        R:=EquationOrder(A);
        tO := OverOrders(R);
        tG := GraphOverOrders(R);
        assert #tO eq #Vertices(tG);
        assert <#Edges(tM), #Vertices(tM), #Vertices(tG), #Edges(tG)> eq data[i];
        // we compare with the minimal inclusions found by testing all inclusions among the
        // overorders
        assert edges(tG) eq min_inclusions_naive(tO);
    end for;

    // The graph must not depend on what was computed before. If the overorders of R are already
    // known when GraphOverOrders is called, for example after WKICM or LoadWKICM, then the minimal
    // overorders of the overorders are not stored, or only some of them are: GraphOverOrders has to
    // compute them.
    fs:=[
          x^6+3*x^5+99*x^2+2997, //1 sing prime
          x^4+4*x^3+6*x^2+44*x+121, //2 sing primes
          (x^2+5)*(x^2+20), //4 sing primes, in a product of two number fields
        ];
    data := [<1, 33, 65>, <2, 16, 26>, <4, 24, 52>]; // number of singular primes, of overorders, of
                                                     // minimal inclusions
    for i->f in fs do
        printf ".";
        // GraphOverOrders is the first computation
        R:=EquationOrder(EtaleAlgebra(f));
        G:=GraphOverOrders(R);
        assert <#SingularPrimes(R), #Vertices(G), #Edges(G)> eq data[i];
        // after WKICM. We use a new algebra, so that no attribute is shared with R.
        R1:=EquationOrder(EtaleAlgebra(f));
        _:=WKICM(R1);
        assert edges(GraphOverOrders(R1)) eq edges(G);
        // after LoadWKICM
        R2:=LoadWKICM(PrintWKICM(R1));
        assert edges(GraphOverOrders(R2)) eq edges(G);
    end for;
    SetAssertions(1);
    printf " all good!";

*/
