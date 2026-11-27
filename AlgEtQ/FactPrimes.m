/////////////////////////////////////////////////////
// Copyright 2025.
// Stefano Marseglia, stefano.marseglia89@gmail.com
// https://stmar89.github.io/index.html
//
// Distributed under the terms of the CC-BY 4.0 licence.
// https://creativecommons.org/licenses/by/4.0/
/////////////////////////////////////////////////////



freeze;

declare verbose FactPrimes, 3;

import "Ord.m" : all_in, zb_dM, zb_inv, hnf, crZQ, MatrixAtoQ;
import "Idl.m" : ideal_from_hnf;

declare attributes AlgEtQOrd :
            SingularPrimes,
            PlacesAboveRationalPrime, // The primes of the maximal order O, built once: an
                                      // AssociativeArray keyed by the rational primes q, whose
                                      // value is the list of the triples <i, P_L, P>, one for each
                                      // prime P of O above q: P is P_L, a prime of the maximal
                                      // order of the i-th component, at the i-th component, and the
                                      // unit ideal at the other components. The list has the
                                      // components in order, and for each of them the primes in the
                                      // order of Factorization(q*O_i). Each P is marked prime, with
                                      // its ramification index and its inertia degree. Written only
                                      // by places_above_rational_prime; read by
                                      // factorizationMaximalOrder, PrimesAbove and
                                      // PlacesAboveRationalPrime. A prime of O is always the same
                                      // object, so what is stored on it (ZBasis, Hash, residue
                                      // field, PrimeBelow, minimal overorders) is computed once.
            PrimesAboveRationalPrime, // The primes of S above the rational primes p: an
                                      // AssociativeArray keyed by p, with value the sequence of the
                                      // primes of S above p, each one a single object. Written only
                                      // by primes_of_order_above; read by PrimesAbove and
                                      // PrimeBelow.
            PrimesBelow; // The cache of PrimeBelow: an AssociativeArray keyed by the orders T
                         // containing S, whose value at T is an AssociativeArray keyed by the
                         // primes P of T, with value PrimeBelow(S,P). The primes are grouped by
                         // their order because ideals over different orders cannot be compared.

///# Primes, factorization and completions

///## Primes
/// Let $S$ be an order in an étale algebra $A$. A `prime` of $S$ is a maximal ideal of $S$, that is, an integral fractional ideal $\mathfrak{p}$ such that $S/\mathfrak{p}$ is a field.
/// It is well-known that a prime $\mathfrak{p}$ of $S$ is invertible if and only if $\mathfrak{p}$ does not contain the conductor of $S$. Such primes are also called `regular`. A prime which is non-invertible is also called `singular`.
/// In particular, there are only finitely many singular primes.

intrinsic IsPrime(I::AlgEtQIdl) -> BoolElt
{Returns whether the given integral fractional ideal is a prime of its order of definition.}
    if not assigned I`IsPrime then
        require IsIntegral(I): "the ideal must be integral";
        // Cheap tests first, from N = [S:I]: S/I is a field only if N = p^f and p is in I (the
        // characteristic of S/I is p); if N = p, then S/I has p elements, so it is a field. Only
        // the remaining cases need the primes above I.
        S:=Order(I);
        is_pp,p,f:=IsPrimePower(Integers()!Index(OneIdeal(S),I));
        if not is_pp or not (Algebra(S)!p in I) then
            bool:=false;
        elif f eq 1 then
            bool:=true;
        else
            prim:=PrimesAbove(I);
            bool:=#prim eq 1 and I eq prim[1];
        end if;
        assert2 bool eq (#PrimesAbove(I) eq 1 and I eq PrimesAbove(I)[1]);
        I`IsPrime:=bool;
    end if;
    return I`IsPrime;
end intrinsic;

///ditto
intrinsic IsMaximal(I::AlgEtQIdl) -> BoolElt
{ " } // "
    return IsPrime(I);
end intrinsic;

///ditto
intrinsic IsMaximalIdeal(I::AlgEtQIdl) -> BoolElt
{ " } // "
    return IsPrime(I);
end intrinsic;

primes_above_p_fp:=function(S,p)
// The primes of an order S above a rational prime p, from the F_p-algebra B = S/pS, without the
// maximal order. The primes of S above p (the maximal ideals containing p) are in bijection with
// the maximal ideals of  B = S/pS, of dimension n: P <-> P/pS. Let rad be the nilradical of B.
// Then B/rad is a product of finite fields k_1 x ... x k_m, and the maximal ideals of  B are the
// preimages of the maximal ideals of B/rad, which are the kernels of the m projections.
// We find them with the Frobenius Fr: x -> x^p, which is F_p-linear on B (char p).
    A:=Algebra(S);
    zb:=ZBasis(S);
    n:=#zb;
    F:=GF(p);
    V:=VectorSpace(F,n); // B = S/pS, in the coordinates of zb

    // 1. The multiplication of B. Let BQ be the matrix of the absolute coordinates of zb (row i is
    //    zb[i]) and R_i the matrix of the multiplication by zb[i] in the absolute basis (row j is
    //    the absolute coordinates of zb[i] times the j-th basis element: RepresentationMatrix).
    //    Then BQ*R_i*BQ^-1 is the matrix of the multiplication by zb[i] in the coordinates of zb:
    //    its row j contains the coordinates of zb[j]*zb[i]. It is an integer matrix, since S is a
    //    ring: we reduce it modulo p. BQ^-1 is the inclusion matrix of S (stored).
    dS,MS:=zb_dM(S); // MS = dS*BQ, an integer matrix
    BQ:=(1/dS)*crZQ(MS);
    BQinv:=zb_inv(S);
    mult_by:=[ ChangeRing(ChangeRing(BQ*RepresentationMatrix(z)*BQinv,Integers()),F) : z in zb ];
    // the product in B of x and y (vectors in V): x = sum x[i]*zb[i], so x*y = sum
    // x[i]*(y*mult_by[i])
    mult:=func< x,y | V!Eltseq(&+[ x[i]*(y*mult_by[i]) : i in [1..n] | x[i] ne 0 ]) >;
    // the identity of B, in the coordinates of zb
    one:=V![ F!c : c in Eltseq(MatrixAtoQ([One(A)])*BQinv) ];

    // 2. The Frobenius Fr: x -> x^p as an n x n matrix over F_p (rows: images of the basis
    //    vectors). zb[i]^p = 1*zb[i]*...*zb[i] (p times), i.e. one*mult_by[i]^p: a power of a
    //    matrix, computed by Magma with repeated squaring.
    Frob:=Matrix(F,n,n,&cat[ Eltseq(one*mult_by[i]^p) : i in [1..n] ]);

    // 3. The nilradical rad = ker(Fr^t) with p^t >= n, and the reduced algebra Q = B/rad, with the
    //    quotient map q: V -> Q. Elements of Q are multiplied through lifts to V.
    t:=1;
    while p^t lt n do
        t+:=1;
    end while;
    rad:=Kernel(Frob^t);
    Q,q:=quo<V|rad>;
    r:=Dimension(Q);
    lift:=func< y | y@@q >;
    multQ:=func< x,y | q(mult(lift(x),lift(y))) >;
    oneQ:=q(one);

    // 4. The fixed points K of the Frobenius on Q (it induces a map on Q, since Fr(rad) c rad):
    //    K = F_p^m, with m the number of primes of S above p.
    FrobQ:=Matrix(F,r,r,&cat[ Eltseq(q(lift(Q.i)*Frob)) : i in [1..r] ]);
    K:=Kernel(FrobQ-IdentityMatrix(F,r));
    m:=Dimension(K);

    // 5. The primitive idempotents of Q: we split the identity with the elements k of a basis of K
    //    (see above). For each current idempotent e we compute the minimal polynomial of k*e in the
    //    algebra eQ, from the first linear dependence among e, k*e, (k*e)^2, ..., and split e
    //    according to its roots (all in F_p, and distinct).
    PX<X>:=PolynomialRing(F);
    idems:=[ oneQ ];
    for k in Basis(K) do
        if #idems eq m then
            break; // all the factors are separated
        end if;
        new:=[];
        for e in idems do
            ke:=multQ(Q!Eltseq(k),e);
            // the powers e, ke, ke^2, ... until they become linearly dependent
            pows:=[ e ];
            repeat
                Append(~pows,multQ(pows[#pows],ke));
            until Rank(Matrix(F,#pows,r,&cat[ Eltseq(v) : v in pows ])) lt #pows;
            // ke^s = c[1]*e + c[2]*ke + ... + c[s]*ke^(s-1): the minimal polynomial of ke in eQ
            c:=Eltseq(Solution(Matrix(F,#pows-1,r,&cat[ Eltseq(v) : v in pows[1..#pows-1] ]),Vector(F,Eltseq(pows[#pows]))));
            minpol:=X^(#c) - &+[ c[i]*X^(i-1) : i in [1..#c] ];
            roots:=[ rt[1] : rt in Roots(minpol) ];
            assert #roots eq Degree(minpol); // distinct roots in F_p, since ke is in F_p^(m_e)
            if #roots eq 1 then
                Append(~new,e); // k is constant on the factors of e: e is not split by k
                continue;
            end if;
            for a in roots do
                // the idempotent prod_{b ne a} (ke - b*e)/(a - b): e on the factors where ke = a
                ea:=e;
                for b in roots do
                    if b ne a then
                        ea:=multQ(ea,ke-b*e)/(a-b);
                    end if;
                end for;
                Append(~new,ea);
            end for;
        end for;
        idems:=new;
    end for;
    assert #idems eq m;

    // 6. The primes. For each primitive idempotent e, the maximal ideal (1-e)*Q of Q, its
    //    preimage U in V (spanned by lifts of a basis of (1-e)*Q and by rad), and the prime of S:
    //    the lattice spanned by p*zb and the lifts of U, i.e. the rows of p*Id and of U in the
    //    coordinates of zb. Its HNF Hc, times MS (= dS times the absolute coordinates of zb), gives
    //    dS times the absolute coordinates of a Z-basis of the prime: one more HNF, and the ideal
    //    is created from it with its hash (ideal_from_hnf in Idl.m), without any product of
    //    elements.
    primes:=[];
    for e in idems do
        maxQ:=sub<Q|[ multQ(Q.i,oneQ-e) : i in [1..r] ]>;
        U:=[ lift(u) : u in Basis(maxQ) ] cat Basis(rad);
        L:=VerticalJoin(p*IdentityMatrix(Integers(),n),
                        Matrix(Integers(),#U,n,&cat[ [ Integers()!c : c in Eltseq(u) ] : u in U ]));
        Hc:=hnf(L); // the prime in the coordinates of zb
        P:=ideal_from_hnf(S,hnf(Hc*MS),dS);
        P`IsPrime:=true;
        Append(~primes,P);
    end for;
    return primes;
end function;

places_above_rational_prime:=function(O, q)
// Given the maximal order O and a rational prime q, returns the list of the triples <i,P_L,P> of
// the primes of O above q (see the attribute PlacesAboveRationalPrime), computing and storing it
// the first time.
    if not assigned O`PlacesAboveRationalPrime then
        O`PlacesAboveRationalPrime:=AssociativeArray();
    end if;
    if not IsDefined(O`PlacesAboveRationalPrime,q) then
        test,OasProd:=IsProductOfOrders(O);
        assert test; // since O is maximal
        tup_one_ideals:=< 1*OL : OL in OasProd >;
        out:=[* *]; // a list: the primes P_L of different components are in different universes
        for i in [1..#OasProd] do
            for pe in Factorization(q*OasProd[i]) do
                tup_p:=tup_one_ideals;
                tup_p[i]:=pe[1]; // we replace the i-th ideal with the prime
                P:=Ideal(O,tup_p);
                ZBasisLLL(P);
                P`IsPrime:=true;
                P`RamificationIndex:=pe[2]; // the exponent of the prime in q*O_i
                P`InertiaDegree:=Ilog(q,Integers()!Index(O,P)); // [O:P] = q^f
                assert2 q^P`InertiaDegree eq Index(O,P);
                Append(~out,<i,pe[1],P>);
            end for;
        end for;
        O`PlacesAboveRationalPrime[q]:=out;
    end if;
    return O`PlacesAboveRationalPrime[q];
end function;

primes_of_order_above:=function(S,p)
// The primes of the order S above the rational prime p from the F_p-algebra S/pS
// (primes_above_p_fp), without the maximal order: used by PrimesAbove when the maximal order of the
// algebra is not known (otherwise the stored primes of the maximal order are contracted, which is
// faster when they are shared by many orders). Computed once and stored in
// S`PrimesAboveRationalPrime. A prime already found by PrimeBelow (an intersection) is returned as
// that object, so that every prime of S is one ideal, whose attributes (ZBasis, hash, residue
// field, minimal overorders) are computed once.
    if not assigned S`PrimesAboveRationalPrime then
        S`PrimesAboveRationalPrime:=AssociativeArray();
    end if;
    if not IsDefined(S`PrimesAboveRationalPrime,p) then
        known:=[];
        if assigned S`PrimesBelow then
            known:=&cat[ [ S`PrimesBelow[T][P] : P in Keys(S`PrimesBelow[T]) ] : T in Keys(S`PrimesBelow) ];
        end if;
        primes:=primes_above_p_fp(S,p);
        for k in [1..#primes] do
            old:=[ Q : Q in known | Q eq primes[k] ];
            if #old gt 0 then
                primes[k]:=old[1];
            else
                ZBasisLLL(primes[k]);
            end if;
        end for;
        S`PrimesAboveRationalPrime[p]:=primes;
    end if;
    return S`PrimesAboveRationalPrime[p];
end function;

/// Given an order $S$ and a prime $P$ of an order $T$ containing $S$, returns the prime $S \cap P$ of $S$ below $P$. The output is stored in $S$: asking again for the same prime returns the same ideal, with all its attributes.
intrinsic PrimeBelow(S::AlgEtQOrd, P::AlgEtQIdl) -> AlgEtQIdl
{Given an order S and a prime P of an order T containing S, returns the prime S meet P of S below P. The output is stored in S: asking again for the same prime returns the same ideal, with all its attributes.}
    T:=Order(P);
    if T eq S then
        return P;
    end if;
    if assigned S`PrimesBelow and IsDefined(S`PrimesBelow,T) and IsDefined(S`PrimesBelow[T],P) then
        return S`PrimesBelow[T][P];
    end if;
    require IsPrime(P) : "the ideal must be prime";
    // The inclusion test is skipped when T is the maximal order, which contains every order: this
    // is the main use, and with many orders S the test would cost more than the intersection
    // itself.
    require (assigned T`IsMaximal and T`IsMaximal) or S subset T : "the order must be contained in the order of the prime";
    // S/(S meet P) embeds in the finite field T/P, so it is a finite domain, hence a field: S meet
    // P is a prime of S, above the rational prime p below P.
    if not assigned S`PrimesBelow then
        S`PrimesBelow:=AssociativeArray();
    end if;
    _,p:=IsPrimePower(Integers()!Index(T,P));
    if assigned S`PrimesAboveRationalPrime and IsDefined(S`PrimesAboveRationalPrime,p) then
        // The primes of S above p are known (primes_of_order_above): S meet P is the one contained
        // in P, since a prime Q of S contained in P is contained in S meet P, hence equal to it (Q
        // is maximal). Q c P is tested on the generators of Q (P is a T-ideal and S c T).
        cands:=[ Q : Q in S`PrimesAboveRationalPrime[p] | all_in(Generators(Q),P) ];
        assert #cands eq 1;
        Q:=cands[1];
        assert2 Q eq OneIdeal(S) meet S!!P;
    else
        Q:=OneIdeal(S) meet S!!P;
        Q`IsPrime:=true;
        // Distinct primes, also of different orders, can be above the same prime of S: if it is
        // already stored, we return it, so that each prime of S is one ideal, whose attributes
        // (residue field, ...) are computed once.
        for T2 in Keys(S`PrimesBelow) do
            for P2 in Keys(S`PrimesBelow[T2]) do
                if S`PrimesBelow[T2][P2] eq Q then
                    Q:=S`PrimesBelow[T2][P2];
                    break T2;
                end if;
            end for;
        end for;
    end if;
    if not IsDefined(S`PrimesBelow,T) then
        S`PrimesBelow[T]:=AssociativeArray();
    end if;
    S`PrimesBelow[T][P]:=Q;
    return Q;
end intrinsic;

/// Given an integral fractional $S$-ideal $I$, returns the sequence of primes of $S$ containing $I$.
intrinsic PrimesAbove(I::AlgEtQIdl) -> SeqEnum[AlgEtQIdl]
{Given an integral fractional S-ideal I, returns the sequence of primes of S containing I.}
    if not assigned I`PrimesAbove then
        require IsIntegral(I): "the ideal must be integral";
        S:=Order(I);
        if assigned I`Factorization then
            primes:=[P[1] : P in I`Factorization];
        elif I eq OneIdeal(S) then
            primes:=[];
        else
            // The primes of S containing I contain N:=MinimalInteger(I), so they are above the
            // rational primes q dividing N; I c P is tested with a few memberships. For each q:
            // - if the primes of S above q are stored (see primes_of_order_above), we use them;
            // - if the maximal order O of the algebra is known, each prime of S above q is
            //   PrimeBelow(S,PO) for a prime PO of O above q, and I c PO meet S if and only if
            //   I c PO. The primes of O above q are built once (places_above_rational_prime) and
            //   shared by all the orders of the algebra: with many orders (e.g. all the overorders)
            //   this is faster than the F_p-algebras below;
            // - otherwise they come from the F_p-algebra S/qS (primes_of_order_above), and the
            //   maximal order is not computed.
            // The output is in a fixed order: by q, then by component of O (or in the order of
            // primes_of_order_above); repetitions (several primes of O above the same prime of S)
            // are removed.
            gI:=Generators(I);
            A:=Algebra(S);
            primes:={@ @};
            for q in PrimeDivisors(MinimalInteger(I)) do
                if assigned S`PrimesAboveRationalPrime and IsDefined(S`PrimesAboveRationalPrime,q) then
                    primes join:={@ P : P in S`PrimesAboveRationalPrime[q] | all_in(gI,P) @};
                elif assigned A`MaximalOrder then
                    O:=MaximalOrder(A);
                    maxS:=IsIdentical(S,O) or (assigned S`IsMaximal and S`IsMaximal);
                    for t in places_above_rational_prime(O,q) do
                        PO:=t[3];
                        // the generators of I over S suffice: S is in O, the order of PO
                        if all_in(gI,PO) then
                            Include(~primes, maxS select PO else PrimeBelow(S,PO));
                        end if;
                    end for;
                else
                    primes join:={@ P : P in primes_of_order_above(S,q) | all_in(gI,P) @};
                end if;
            end for;
            primes:=IndexedSetToSequence(primes);
            assert2 forall{P : P in primes | IsIntegral(P)};
            assert2 forall{P : P in primes | I subset P};
            assert2 forall{P : P in primes | IsPrimePower(Integers() ! Index(S,P))};
        end if;
        for P in primes do
            ZBasisLLL(P);
        end for;
        I`PrimesAbove:=primes;
    end if;
    return I`PrimesAbove;
end intrinsic;

intrinsic PlacesAboveRationalPrime(E::AlgEtQ,p::RngIntElt)->SeqEnum[AlgEtQIdl]
{Given an étale algebra and a rational prime, returns the primes of maximal order of the algebra containing the rational prime.}
    require IsPrime(p) : "The integer p needs to be a prime number";
    // stored in the maximal order, see the attribute PlacesAboveRationalPrime
    return [ t[3] : t in places_above_rational_prime(MaximalOrder(E),p) ];
end intrinsic;

intrinsic SingularPrimes(R::AlgEtQOrd) -> SeqEnum[AlgEtQIdl]
{Returns the non-invertible primes of the order.}
    if not assigned R`SingularPrimes then
        // The conductor is computed with the maximal order, before PrimesAbove: so PrimesAbove
        // contracts the stored primes of the maximal order, which is faster than the F_p-algebras
        // when the singular primes of many orders are needed (see PrimesAbove).
        R`SingularPrimes:=PrimesAbove(Conductor(R));
    end if;
    return R`SingularPrimes;
end intrinsic;

///ditto
intrinsic NonInvertiblePrimes(R::AlgEtQOrd) -> SeqEnum[AlgEtQIdl]
{ " } // "
    return SingularPrimes(R);
end intrinsic;

///## Factorization
/// Every integral fractional $S$-ideal which is coprime to the conductor of $S$ can be written as a product of invertible primes of $S$.
/// This factorization is unique up to permutation of the factors.
///
/// Consider the decomposition $K_1\times \cdots \times K_n$ of $A$ into the direct product of its components.
/// Then every prime $\mathfrak{P}$ of the maximal order $\mathcal{O}_A$ of $A$ is generated by a maximal ideal in exactly one component and by $1$ in every other component.
/// It follows that the factorization problem for fractional ideals over the maximal order can be immediately reduced to the number field case.
/// We deal the case of non-maximal orders by taking intersections.

function factorizationMaximalOrder(I)
//given an ideal of the maximal order of an étale algebra, returns the factorization into a product
//of prime ideals
    O:=Order(I);
    assert IsMaximal(O);
    test,IasProd:=IsProductOfIdeals(I);
    assert test;
    fac:=[]; //this will be the factorization of I
    for i in [1..#IasProd] do
        IL:=IasProd[i];
        if not One(Order(IL)) in IL then
            facL:=Factorization(IL); // < <P,e> : ... >;
            assert #facL gt 0;
            for p in facL do
                // the prime of O which is p[1] at the i-th component, built once (see
                // places_above_rational_prime)
                Ps:=[ t[3] : t in places_above_rational_prime(O,Minimum(p[1])) | t[1] eq i and t[2] eq p[1] ];
                assert #Ps eq 1;
                Append(~fac,<Ps[1],p[2]>);
            end for;
         end if;
    end for;
    // The powers are formed by repeated products, not with '^': for an invertible ideal '^' calls
    // TwoGeneratingSet, which can call CoprimeRepresentative, PrimesAbove and Factorization of an
    // equal ideal again, an infinite recursion at assertion level 2 (Magma then crashes).
    assert2 I eq &*[ &*[ p[1] : i in [1..p[2]] ] : p in fac ];
    return fac;
end function;

/// Given an integral fractional $S$-ideal $I$ coprime with the conductor of $S$ (hence invertible in $S$), returns its factorization into a product of primes of $S$.
intrinsic Factorization(I::AlgEtQIdl) -> Tup
{Given an integral fractional S-ideal I coprime with the conductor of S (hence invertible in S), returns its factorization into a product of primes of S.}
    if not assigned I`Factorization then
        S:=Order(I);
        require IsIntegral(I) and I ne OneIdeal(S): "the argument must be a proper integral ideal";
        if IsMaximal(S) then
            I`Factorization:=factorizationMaximalOrder(I);
        else
            fS:=Conductor(S);
            require IsCoprime(fS,I): "the ideal must be coprime with the conductor of the order of definition";
            O:=MaximalOrder(Algebra(I));
            IO:=O !! I;
            facO:=factorizationMaximalOrder(IO);
            // I is coprime with the conductor of S. So every prime P of S containing I is below
            // exactly one prime PO of O containing IO, and the valuations of I at P and of IO at PO
            // are equal, since S and O are equal locally at P. Hence we only need to intersect each
            // PO with S.
            facS:=<>;
            for pO in facO do
                P:=PrimeBelow(S,pO[1]);
                ZBasisLLL(P);
                P`IsPrime:=true;
                Append(~facS, <P,pO[2]>);
            end for;
            assert2 #{ p[1] : p in facS } eq #facS; // the primes of S are distinct
            // no '^': see factorizationMaximalOrder
            assert2 I eq &*[ &*[ p[1] : i in [1..p[2]] ] : p in facS ];
            I`Factorization:=facS;
        end if;
     end if;
     return I`Factorization;
end intrinsic;


/* TESTS

    printf "### Testing Primes and Factorizaton:";
    //AttachSpec("../spec");
    SetAssertions(2);

    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16);
    A:=EtaleAlgebra(f);
    E1:=EquationOrder(A);
    ff:=Conductor(E1);
    assert PrimesAbove(Conductor(E1)) eq SingularPrimes(E1);
    printf ".";

    f:=x^6 + 8*x^5 + 50*x^4 + 200*x^3 + 1250*x^2 + 5000*x + 15625;
    A:=EtaleAlgebra(f);
    assert #SingularPrimes(MaximalOrder(A)) eq 0;
    R:=EquationOrder(A);
    assert #SingularPrimes(R) eq 5;
    assert #NonInvertiblePrimes(R) eq 5;
    printf ".";

    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    E1:=EquationOrder(A);
    E2:=ProductOfEquationOrders(A);

    _:=PrimesAbove(Conductor(E1));
    _:=PrimesAbove(Conductor(E2));
    assert IsGorenstein(E1);
    assert IsGorenstein(E2);

    // random principal ideals: an ideal generated by many random elements is almost always E1
    ids:=[ Random(E1)*E1 : i in [1..5] ];
    ids0:=[ I : I in ids | I ne OneIdeal(E1) and IsInvertible(I) ];
    ids:=[];
    for I in ids0 do
        _,J:=CoprimeRepresentative(I,Conductor(E1));
        // J can be E1 itself, for example when the inverse of a generator of I is short
        if J ne OneIdeal(E1) then
            Append(~ids,J);
        end if;
    end for;
    facs:=[* Factorization(I) : I in ids *]; // a list: tuples of different lengths are not in the
                                             // same universe
    // an ideal with many prime factors, coprime with the conductor.
    // We compare its primes with PrimesAbove of a copy without the attribute Factorization.
    ps:=[ p : p in PrimesUpTo(20) | IsCoprime(p*OneIdeal(E1),Conductor(E1)) ];
    I:=(&*ps)*OneIdeal(E1);
    fac:=Factorization(I);
    J:=Ideal(E1,ZBasis(I));
    assert Seqset(PrimesAbove(J)) eq { p[1] : p in fac };
    printf ".";

    // testing PlacesAboveRationalPrime: the primes above different rational primes are stored
    // separately
    _<x>:=PolynomialRing(Integers());
    A:=EtaleAlgebra((x^2+5)*(x^2+7));
    O:=MaximalOrder(A);
    P2:=PlacesAboveRationalPrime(A,2);
    P3:=PlacesAboveRationalPrime(A,3);
    assert Seqset(P2) eq Seqset(PrimesAbove(2*O));
    assert Seqset(P3) eq Seqset(PrimesAbove(3*O));
    assert Seqset(P2) ne Seqset(P3);
    assert PlacesAboveRationalPrime(A,2) eq P2; // the stored value
    assert Keys(O`PlacesAboveRationalPrime) eq {2,3};
    // the primes of O are built once: the same objects in PlacesAboveRationalPrime, PrimesAbove and
    // Factorization, marked prime, with their ramification indices and inertia degrees
    assert forall{ P : P in PrimesAbove(2*O) cat [ f[1] : f in Factorization(6*O) ] |
                        exists{ Q : Q in P2 cat P3 | IsIdentical(P,Q) } };
    for pP in [ <2,P2>, <3,P3> ] do
        p,Ps:=Explode(pP);
        assert forall{ P : P in Ps | assigned P`IsPrime and assigned P`RamificationIndex and assigned P`InertiaDegree };
        assert forall{ P : P in Ps | P`RamificationIndex eq Valuation(A!p,P) and p^P`InertiaDegree eq Index(O,P) };
        assert &+[ RamificationIndex(P)*InertiaDegree(P) : P in Ps ] eq AbsoluteDimension(A);
    end for;
    printf ".";

    // PrimeBelow: the prime S meet P of S below a prime P of an order T containing S. It is stored:
    // the same ideal object is returned again, also for an equal prime which is a different object.
    // Primes of the maximal order and of an intermediate order T, in an algebra with two
    // components.
    A:=EtaleAlgebra((x^2+5)*(x^3-2000));
    E:=EquationOrder(A);
    O:=MaximalOrder(A);
    T:=[ S : S in OverOrders(E) | S ne E and S ne O ][1];
    for P in PrimesAbove(Conductor(E)) do
        for U in [ O, T ] do
            for Q in PrimesAbove(U!!P) do
                Q2:=PrimeBelow(E,Q);
                assert Order(Q2) eq E and IsPrime(Q2) and Q2 eq OneIdeal(E) meet E!!Q;
                assert Q2 eq P; // E has a unique prime below Q, and P is a prime of E below Q
                assert IsIdentical(PrimeBelow(E,Q),Q2);
                Q3:=Ideal(U,ZBasis(Q)); // equal to Q, a different object
                assert not IsIdentical(Q,Q3) and IsIdentical(PrimeBelow(E,Q3),Q2);
            end for;
        end for;
    end for;
    // one ideal for each prime of E: the primes of E below the different primes of O and of T above
    // the same prime P of E are all the same object, which is P itself (PrimesAbove uses
    // PrimeBelow)
    for P in PrimesAbove(Conductor(E)) do
        below:=[ PrimeBelow(E,Q) : Q in PrimesAbove(O!!P) cat PrimesAbove(T!!P) ];
        assert forall{ Q2 : Q2 in below | IsIdentical(Q2,P) };
    end for;
    // P of T itself
    assert IsIdentical(PrimeBelow(T,PrimesAbove(Conductor(T))[1]),PrimesAbove(Conductor(T))[1]);
    try
        _:=PrimeBelow(O,PrimesAbove(Conductor(E))[1]); // O is not contained in E
        ok:=false;
    catch err
        ok:=true;
    end try;
    assert ok;
    printf ".";
    // at assertion level 2, the check of a factorization formed the prime powers with '^', which
    // for this ideal led to an infinite recursion (TwoGeneratingSet, CoprimeRepresentative,
    // PrimesAbove, Factorization) and a crash of Magma; found by the timing test PicardGroup_1.8 at
    // level 2
    SetAssertions(2);
    _<x>:=PolynomialRing(Integers());
    A:=EtaleAlgebra(x^6-3*x^5-3*x^4+65*x^3-48*x^2-768*x+4096);
    assert [ DefiningPolynomial(K) : K in Components(A) ] eq [ x^2-7*x+16, x^4+4*x^3+9*x^2+64*x+256 ];
    zbS:=[ <[ 0, 0 ], [ -1/9, 1/48, 1/12, 11/144 ]>, <[ 0, 0 ], [ 0, -1/24, 1/6, -1/24 ]>,
           <[ 0, 0 ], [ 1/9, 3/16, 1/12, 19/144 ]>, <[ 0, 0 ], [ 1/9, -1/6, 0, 5/18 ]>,
           <[ 1, 0 ], [ 0, 0, 0, 0 ]>, <[ 0, 1 ], [ 0, 0, 0, 0 ]> ];
    gI:=[ <[ 2, 0 ], [ 0, 0, 0, 0 ]>, <[ 3, 1 ], [ 0, 0, 0, 0 ]>, <[ 0, 0 ], [ 1, 0, 0, 0 ]> ];
    S:=Order([ A!t : t in zbS ]);
    I:=Ideal(S,[ A!t : t in gI ]);
    TwoGeneratingSet(I);
    assert #Generators(I) le 2 and I eq Ideal(S,Generators(I));
    printf ".";
    // IsPrime decides from the index first: fresh copies, so that no flag is set
    A:=EtaleAlgebra((x^2+5)*(x^2-2));
    E:=EquationOrder(A);
    for P in PrimesAbove(6*OneIdeal(E)) do
        assert IsPrime(Ideal(E,ZBasis(P))); // index p or p^f
    end for;
    assert not IsPrime(Ideal(E,ZBasis(6*OneIdeal(E)))); // the index is not a prime power
    assert not IsPrime(Ideal(E,ZBasis(PrimesAbove(3*OneIdeal(E))[1]^2))); // index p^2, not prime
    assert not IsPrime(Ideal(E,ZBasis(9*OneIdeal(E)))); // 3 is not in 9E
    printf ".";
    // Without the maximal order: PrimesAbove from the F_p-algebras S/pS. The same primes (as
    // lattices) as through the maximal order, and the same objects afterwards, also from
    // PrimeBelow.
    for f in [ (x^2+5)*(x^2+7)*(x^2+11), x^6+3*x^5+99*x^2+2997, x^4+16 ] do
        A:=EtaleAlgebra(f); E:=EquationOrder(A);
        pp:=PrimesAbove(30*OneIdeal(E));
        assert not assigned A`MaximalOrder;
        B:=EtaleAlgebra(f); EB:=EquationOrder(B); _:=MaximalOrder(B);
        assert { myHash(P) : P in pp } eq { myHash(P) : P in PrimesAbove(30*OneIdeal(EB)) };
        O:=MaximalOrder(A);
        pp2:=PrimesAbove(6*OneIdeal(E)); // now the maximal order is known
        assert forall{ P : P in pp2 | exists{ Q : Q in pp | IsIdentical(P,Q) } };
        assert forall{ PO : PO in PrimesAbove(30*O) | exists{ Q : Q in pp | IsIdentical(PrimeBelow(E,PO),Q) } };
    end for;
    printf ".";
    SetAssertions(1);
    printf " all good!";

*/
