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
            MultiplicationMatrices, // The matrices of the multiplications by the elements of
                                    // ZBasis(S) in the coordinates of ZBasis(S) (integer matrices,
                                    // for the working Z-basis they were computed with: a tuple
                                    // <matrices, ZBasis(S)>). Written and read only by
                                    // fp_multiplication.
            PrimesBelow; // The cache of PrimeBelow: an AssociativeArray keyed by the orders T
                         // containing S, whose value at T is an AssociativeArray keyed by the
                         // primes P of T, with value PrimeBelow(S,P). The primes are grouped by
                         // their order because ideals over different orders cannot be compared.

fp_multiplication:=function(S,p)
// The F_p-algebra S/pS in the coordinates of zb = ZBasis(S): returns F = GF(p), V = F^n, the
// matrices mult_by[i] of the multiplications by zb[i] (row j: the coordinates of zb[j]*zb[i]), and
// the identity one, as a vector of V. Let BQ be the matrix of the absolute coordinates of zb and
// R_i the matrix of the multiplication by zb[i] in the absolute basis (RepresentationMatrix): then
// BQ*R_i*BQ^-1 is the matrix of the multiplication by zb[i] in the coordinates of zb, an integer
// matrix since S is a ring, which we reduce modulo p. BQ^-1 is the inclusion matrix of S.
// The integer matrices do not depend on p: they are computed once and stored in
// S`MultiplicationMatrices, together with the Z-basis they refer to (if the working Z-basis is
// replaced, they are computed again).
    A:=Algebra(S);
    zb:=ZBasis(S);
    n:=#zb;
    F:=GF(p);
    V:=VectorSpace(F,n);
    BQinv:=zb_inv(S);
    if not assigned S`MultiplicationMatrices or not IsIdentical(S`MultiplicationMatrices[2],zb) then
        dS,MS:=zb_dM(S); // MS = dS*BQ, an integer matrix
        BQ:=(1/dS)*crZQ(MS);
        S`MultiplicationMatrices:=<[ ChangeRing(BQ*RepresentationMatrix(z)*BQinv,Integers()) : z in zb ],zb>;
    end if;
    mult_by:=[ ChangeRing(M,F) : M in S`MultiplicationMatrices[1] ];
    one:=V![ F!c : c in Eltseq(MatrixAtoQ([One(A)])*BQinv) ];
    return F,V,mult_by,one;
end function;

is_field_fp:=function(I,p)
// Given an integral ideal I of an order S with p in I (so that B = S/I is an F_p-algebra), returns
// whether B is a field. A finite commutative F_p-algebra is reduced if and only if the Frobenius x
// -> x^p (an F_p-linear map) is injective, since x^p = 0 for x ne 0 exactly when B has nilpotents;
// a reduced B is a product of finite fields, and the fixed points of the Frobenius are F_p in each
// factor: B is a field if and only if they form a space of dimension 1. B = V/U with V = S/pS (see
// fp_multiplication) and U = I/pS, spanned by the coordinates (integers) of ZBasis(I) modulo p. The
// Frobenius of B comes from the lifts x of the basis of B: x^p = one*L_x^p, with L_x the matrix of
// the multiplication by x, a sum of the mult_by[i].
    S:=Order(I);
    F,V,mult_by,one:=fp_multiplication(S,p);
    n:=Dimension(V);
    U:=sub<V| [ V![ F!(Integers()!c) : c in r ] : r in AbsoluteCoordinates(ZBasis(I),S) ]>;
    B,q:=quo<V|U>;
    f:=Dimension(B);
    images:=[];
    for i in [1..f] do
        x:=B.i@@q;
        Lx:=&+[ x[k]*mult_by[k] : k in [1..n] | x[k] ne 0 ];
        Append(~images,q(one*Lx^p));
    end for;
    FrobB:=Matrix(F,f,f,&cat[ Eltseq(y) : y in images ]);
    return Dimension(Kernel(FrobB)) eq 0 and Dimension(Kernel(FrobB-IdentityMatrix(F,f))) eq 1;
end function;

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
        // the remaining cases need a computation in the F_p-algebra S/I.
        S:=Order(I);
        is_pp,p,f:=IsPrimePower(Integers()!Index(OneIdeal(S),I));
        if assigned I`Factorization then
            // a known factorization decides (see also Factorization)
            bool:=#I`Factorization eq 1 and I`Factorization[1][2] eq 1;
        elif not is_pp or not (Algebra(S)!p in I) then
            bool:=false;
        elif f eq 1 then
            bool:=true;
        else
            // S/I is an F_p-algebra of dimension f >= 2: a direct test, without computing the
            // primes (see is_field_fp)
            bool:=is_field_fp(I,p);
        end if;
        assert2 bool eq (#PrimesAbove(I) eq 1 and I eq PrimesAbove(I)[1]);
        I`IsPrime:=bool;
        if bool and not assigned I`MinimalInteger then
            I`MinimalInteger:=p; // I meet Z = pZ
        end if;
        // A prime which is invertible (e.g. a prime of the maximal order) is coprime to the
        // conductor: its factorization is itself, and so are the primes above it.
        if bool and ((assigned S`IsMaximal and S`IsMaximal) or (assigned I`IsInvertible and I`IsInvertible)) then
            if not assigned I`Factorization then
                I`Factorization:=< <I,1> >;
            end if;
            if not assigned I`PrimesAbove then
                I`PrimesAbove:=[ I ];
            end if;
        end if;
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
    // 1. The multiplication of B: the matrices mult_by[i] of the multiplications by zb[i] in the
    //    coordinates of zb, and the identity one (see fp_multiplication).
    F,V,mult_by,one:=fp_multiplication(S,p); // B = S/pS, in the coordinates of zb
    dS,MS:=zb_dM(S); // MS = dS times the matrix of the absolute coordinates of zb
    // the product in B of x and y (vectors in V): x = sum x[i]*zb[i], so x*y = sum
    // x[i]*(y*mult_by[i])
    mult:=func< x,y | V!Eltseq(&+[ x[i]*(y*mult_by[i]) : i in [1..n] | x[i] ne 0 ]) >;

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
    Frobt:=Frob^t;
    rad:=Kernel(Frobt);
    Q,q:=quo<V|rad>;
    r:=Dimension(Q);
    lift:=func< y | y@@q >;
    oneQ:=q(one);

    // 4. The fixed points K of the Frobenius on Q (it induces a map on Q, since Fr(rad) c rad):
    //    K = F_p^m, with m the number of primes of S above p.
    FrobQ:=Matrix(F,r,r,&cat[ Eltseq(q(lift(Q.i)*Frob)) : i in [1..r] ]);
    K:=Kernel(FrobQ-IdentityMatrix(F,r));
    m:=Dimension(K);

    // 5. The decomposition of Q into its m field factors k_j = e_j*Q, with matrices: for an element
    //    k of K, the multiplication by k on Q is diagonalizable with eigenvalues in F_p (k is in
    //    F_p^m), and its eigenspaces are the sums of the factors on which k takes the same value.
    //    Intersecting the eigenspaces of the elements of a basis of K separates all the factors,
    //    since K = F_p^m is spanned by its basis. Then e_j is the component of 1 in the j-th
    //    factor, and the j-th maximal ideal (1-e_j)*Q is the sum of the other factors.
    // the matrix of the multiplication by x in Q (rows: the products of the basis vectors with x)
    mult_by_in_Q:=function(x)
        y:=lift(x);
        M:=&+[ y[i]*mult_by[i] : i in [1..n] | y[i] ne 0 ]; // v -> v*y on V
        return Matrix(F,r,r,&cat[ Eltseq(q(lift(Q.i)*M)) : i in [1..r] ]);
    end function;
    Qr:=VectorSpace(F,r); // Q in its coordinates
    blocks:=[ Qr ];
    for k in Basis(K) do
        if #blocks eq m then
            break; // all the factors are separated
        end if;
        Lk:=mult_by_in_Q(Q!Eltseq(k));
        roots:=[ rt[1] : rt in Roots(MinimalPolynomial(Lk)) ];
        if #roots gt 1 then
            eig:=[ Kernel(Lk - a*IdentityMatrix(F,r)) : a in roots ];
            blocks:=[ W meet E : W in blocks, E in eig | Dimension(W meet E) gt 0 ];
        end if;
    end for;
    assert #blocks eq m;
    // the components of 1 in the factors, and the maximal ideals
    big:=Matrix(F,r,r,&cat[ &cat[ Eltseq(b) : b in Basis(W) ] : W in blocks ]); // the bases, joined
    c:=Eltseq(Solution(big,Vector(F,Eltseq(oneQ))));
    starts:=[ 1 + &+[ Integers() | Dimension(blocks[i]) : i in [1..j-1] ] : j in [1..m] ];
    idems:=[ Q!Eltseq(&+[ c[starts[j]+i-1]*big[starts[j]+i-1] : i in [1..Dimension(blocks[j])] ]) : j in [1..m] ];
    maxQs:=[* *]; // a list: subspaces of different dimensions
    for j in [1..m] do
        gensM:=[ Q | ];
        for i in [1..m] do
            if i ne j then
                gensM cat:=[ Q!Eltseq(b) : b in Basis(blocks[i]) ];
            end if;
        end for;
        Append(~maxQs,sub<Q|gensM>);
    end for;

    // 6. The primes. For each primitive idempotent e, the maximal ideal (1-e)*Q of Q, its
    //    preimage U in V (spanned by lifts of a basis of (1-e)*Q and by rad), and the prime of S:
    //    the lattice spanned by p*zb and the lifts of U, i.e. the rows of p*Id and of U in the
    //    coordinates of zb. Its HNF Hc, times MS (= dS times the absolute coordinates of zb), gives
    //    dS times the absolute coordinates of a Z-basis of the prime: one more HNF, and the ideal
    //    is created from it with its hash (ideal_from_hnf in Idl.m), without any product of
    //    elements.
    // 7. Few generators for each prime, instead of the n elements of its Z-basis: the operations
    //    on an ideal, e.g. the colon (P:P) (one multiplication matrix per generator) or P*P, cost
    //    in proportion to the number of its generators. B is the product of the local algebras
    //    e_j*B, for the idempotents e_j of B lifting the primitive idempotents of Q: a lift x of
    //    e_j(mod rad) gives e_j = x^(p^t) = x*Fr^t, since (e+r)^p = e+r^p in characteristic p and
    //    r^(p^t) = 0 for r in rad. The maximal ideal of B corresponding to e_j is
    //    (1-e_j)*B + m_j, with m_j = e_j*rad the maximal ideal of e_j*B. If v_1,...,v_g generate
    //    m_j as a B-module, then (1-e_j)+v_1, v_2, ..., v_g generate (1-e_j)*B + m_j (since
    //    B*((1-e_j)+v) = (1-e_j)*B + B*v, the two idempotents being orthogonal), and 1-e_j alone
    //    if m_j = 0. The prime is generated over S by p and lifts of these elements. The v_i are
    //    chosen greedily, each one outside the B-module spanned by the previous ones (B*v is the
    //    row space of the matrix of the multiplication by v); g is small: it is the number of
    //    generators of the maximal ideal of S_P/pS_P, 1 or 0 if P is invertible.
    // rows: zb[i]*v
    mult_matrix:=func< v | Matrix(F,n,n,&cat[ Eltseq(v*mult_by[i]) : i in [1..n] ]) >;
    rad_mat:=Matrix(F,Dimension(rad),n,&cat[ Eltseq(b) : b in Basis(rad) ]);
    to_S:=func< vs | [ A | &+[ (Integers()!c[i])*zb[i] : i in [1..n] | c[i] ne 0 ] : c in vs ] >;
    primes:=[];
    for j in [1..m] do
        e:=idems[j];
        maxQ:=maxQs[j];
        U:=[ lift(u) : u in Basis(maxQ) ] cat Basis(rad);
        L:=VerticalJoin(p*IdentityMatrix(Integers(),n),
                        Matrix(Integers(),#U,n,&cat[ [ Integers()!c : c in Eltseq(u) ] : u in U ]));
        Hc:=hnf(L); // the prime in the coordinates of zb
        P:=ideal_from_hnf(S,hnf(Hc*MS),dS);
        P`IsPrime:=true;
        P`MinimalInteger:=p; // P meet Z = pZ
        // the generators, see 7.
        eB:=lift(e)*Frobt; // the idempotent of B lifting e
        assert mult(eB,eB) eq eB;
        mj:=sub<V| Dimension(rad) gt 0 select Rows(rad_mat*mult_matrix(eB)) else [] >; // e*rad
        if Dimension(mj) eq 0 then
            // The local factor e*B of B is a field: P*S_P = p*S_P is principal, so P is invertible,
            // with multiplicator ring S. This is free here, while the test (P:P) eq S costs a colon
            // ideal and the creation of an order: e.g. above the primes dividing the resultants of
            // the components (which divide the discriminant with exponent >= 2) most primes are of
            // this kind, and SingularPrimes tests them all. (m_j ne 0 does not decide: Z[sqrt(-3)]
            // is singular at the prime above 2, Z[i] is not.)
            P`IsInvertible:=true;
            P`MultiplicatorRing:=S;
            assert2 MultiplicatorRing(Ideal(S,ZBasis(P))) eq S;
        end if;
        vs:=[];
        W:=sub<V|>;
        for v in Basis(mj) do
            if not v in W then
                Append(~vs,v);
                W:=W+sub<V|Rows(mult_matrix(v))>; // B*v, contained in mj
            end if;
            if W eq mj then
                break;
            end if;
        end for;
        gensB:=#vs eq 0 select [ one-eB ] else [ one-eB+vs[1] ] cat vs[2..#vs];
        // If p is inert in S (B is a field: e = 1 and m_j = 0), then 1-e = 0 and P = pS is
        // generated by p alone.
        gensB:=[ g : g in gensB | not IsZero(g) ];
        gens:=[ A!p ] cat to_S(gensB);
        assert2 Ideal(S,gens) eq P;
        P`Generators:=gens;
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
                // if this prime of O was found from O/qO (see primes_of_order_above, e.g. when O
                // was found maximal by IsMaximal), we use that object: each prime is one ideal
                if assigned O`PrimesAboveRationalPrime and IsDefined(O`PrimesAboveRationalPrime,q) then
                    old:=[ Q : Q in O`PrimesAboveRationalPrime[q] | Q eq P ];
                    if #old gt 0 then
                        P:=old[1];
                    end if;
                end if;
                ZBasisLLL(P);
                P`IsPrime:=true;
                P`MinimalInteger:=q; // P meet Z = qZ
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
// (primes_above_p_fp), without the maximal order: used by SingularPrimes, and by PrimesAbove when
// the maximal order of the algebra is not known (otherwise the stored primes of the maximal order
// are contracted, which is faster when they are shared by many orders). Computed once and stored in
// S`PrimesAboveRationalPrime. A prime already found by PrimeBelow (an intersection), or by
// places_above_rational_prime for a maximal order, is returned as that object, so that every prime
// of S is one ideal, whose attributes (ZBasis, hash, residue field, minimal overorders) are
// computed once.
    if not assigned S`PrimesAboveRationalPrime then
        S`PrimesAboveRationalPrime:=AssociativeArray();
    end if;
    if not IsDefined(S`PrimesAboveRationalPrime,p) then
        known:=[];
        if assigned S`PrimesBelow then
            known:=&cat[ [ S`PrimesBelow[T][P] : P in Keys(S`PrimesBelow[T]) ] : T in Keys(S`PrimesBelow) ];
        end if;
        // the primes of a maximal order S built from the factorization of p in the components (see
        // places_above_rational_prime)
        if assigned S`PlacesAboveRationalPrime and IsDefined(S`PlacesAboveRationalPrime,p) then
            known cat:=[ t[3] : t in S`PlacesAboveRationalPrime[p] ];
        end if;
        primes:=primes_above_p_fp(S,p);
        // Their Z-bases are LLL-reduced only when they are returned (by PrimesAbove and
        // SingularPrimes): many of them are only tested here.
        for k in [1..#primes] do
            old:=[ Q : Q in known | Q eq primes[k] ];
            if #old gt 0 then
                primes[k]:=old[1];
            end if;
        end for;
        S`PrimesAboveRationalPrime[p]:=primes;
    end if;
    return S`PrimesAboveRationalPrime[p];
end function;

/// Given an order $S$ and a rational prime $p$, returns the primes of $S$ above $p$, that is, the maximal ideals of $S$ containing $p$. If $S$ is known to be maximal, they come from the factorization of $p$ in the components of the algebra, with their ramification indices and inertia degrees; otherwise from the $\mathbb{F}_p$-algebra $S/pS$, without computing the maximal order. They are computed once and stored in $S$.
intrinsic PrimesAboveRationalPrime(S::AlgEtQOrd, p::RngIntElt) -> SeqEnum[AlgEtQIdl]
{Given an order S and a rational prime p, returns the primes of S above p, that is, the maximal ideals of S containing p. If S is known to be maximal, they come from the factorization of p in the components of the algebra, with their ramification indices and inertia degrees; otherwise from the F_p-algebra S/pS, without computing the maximal order. They are computed once and stored in S.}
    require IsPrime(p) : "p must be a prime number";
    if assigned S`IsMaximal and S`IsMaximal then
        return [ t[3] : t in places_above_rational_prime(S,p) ];
    end if;
    return primes_of_order_above(S,p);
end intrinsic;

store_as_maximal:=procedure(S)
// S is found to be maximal (it has no singular primes): it is stored as such and, unless the
// maximal order is already known, as the maximal order of the algebra, so that MaximalOrder returns
// it without computing the maximal orders of the components. The orders of the components, if
// already built (see IsProductOfOrders), are declared maximal to Magma.
    S`IsMaximal:=true;
    S`SingularPrimes:=[];
    A:=Algebra(S);
    if not assigned A`MaximalOrder then
        A`MaximalOrder:=S;
    end if;
    if assigned S`IsProductOfOrders and S`IsProductOfOrders[1] then
        for OL in S`IsProductOfOrders[2] do
            SetOrderMaximal(OL,true);
        end for;
    end if;
end procedure;

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
    p:=MinimalInteger(P); // P meet Z = pZ (stored for the primes built by the package)
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
        Q`MinimalInteger:=p; // Q meet Z = P meet Z = pZ
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
        // Without the maximal order O (and the conductor (R:O)): a prime P of R above p is singular
        // if and only if (P:P) ne R (if P is not invertible then P*P^-1 = P, so P^-1, which is
        // strictly bigger than R, is contained in (P:P)). If P is singular then R is not maximal at
        // P, so p divides [O:R], and p^2 divides disc(R) = [O:R]^2*disc(O). So the candidates are
        // the primes of R above the p with p^2 | disc(R) (from the F_p-algebras R/pR, see
        // primes_of_order_above), and the multiplicator rings (P:P) computed here are needed later
        // anyway (e.g. MinimalOverOrdersAtPrime). The primes whose local factor of R/pR is a field
        // are marked invertible by primes_above_p_fp, so that for them the test costs nothing. The
        // output is in a fixed order: by p, then in the order of primes_of_order_above.
        sing:=[];
        if not (assigned R`IsMaximal and R`IsMaximal) then
            for pe in Factorization(Discriminant(R)) do
                p,e:=Explode(pe);
                if e ge 2 then
                    // the primes of R/pR are marked prime: IsInvertible tests (P:P) eq R, and
                    // the order (P:P) is created only for a singular P
                    sing cat:=[ P : P in primes_of_order_above(R,p) | not IsInvertible(P) ];
                end if;
            end for;
        end if;
        assert2 { myHash(P) : P in sing } eq { myHash(P) : P in PrimesAbove(Conductor(R)) };
        for P in sing do
            ZBasisLLL(P); // as PrimesAbove does for its output
        end for;
        R`SingularPrimes:=sing;
        // the same answer as IsMaximal, whichever is called first
        if #sing eq 0 then
            store_as_maximal(R);
        else
            R`IsMaximal:=false;
        end if;
    end if;
    return R`SingularPrimes;
end intrinsic;

/// Returns whether the given order is the maximal order of the étale algebra. If the maximal order is not known, it is not computed: an order is maximal if and only if it has no singular primes, which are tested from the F_p-algebras $S/pS$ (see `SingularPrimes`). If the order is maximal, it is stored as the maximal order of the algebra.
intrinsic IsMaximal(S::AlgEtQOrd) -> BoolElt
{Returns whether the given order is the maximal order of the étale algebra. If the maximal order is not known, it is not computed: an order is maximal if and only if it has no singular primes, which are tested from the F_p-algebras S/pS (see SingularPrimes). If the order is maximal, it is stored as the maximal order of the algebra.}
    if not assigned S`IsMaximal then
        A:=Algebra(S);
        if assigned A`MaximalOrder then
            S`IsMaximal:=S eq A`MaximalOrder;
        elif assigned S`SingularPrimes then
            S`IsMaximal:=#S`SingularPrimes eq 0;
        else
            // As in SingularPrimes, but stopping at the first singular prime: the candidates are
            // the primes above the p with p^2 | disc(S), and a prime P is singular if and only if
            // it is not invertible, i.e. (P:P) ne S (see IsInvertible). The callers which need the
            // maximal order anyway compute it first (e.g. PicardGroup), and then this is only a
            // comparison.
            maximal:=true;
            for pe in Factorization(Discriminant(S)) do
                if pe[2] ge 2 and exists{ P : P in primes_of_order_above(S,pe[1]) | not IsInvertible(P) } then
                    maximal:=false;
                    break;
                end if;
            end for;
            S`IsMaximal:=maximal;
            if maximal then
                // stored as the maximal order of the algebra (IsProductOfOrders builds the orders
                // of the components from ZBasis(S) when they are needed)
                store_as_maximal(S);
            end if;
        end if;
    end if;
    return S`IsMaximal;
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
        // The maximal order is needed in both cases below (its components, or O itself): computed
        // first, so that IsMaximal(S) is only a comparison (see IsMaximal).
        _:=MaximalOrder(Algebra(S));
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
        // the factorization decides whether I is prime (see also IsPrime)
        if not assigned I`IsPrime then
            I`IsPrime:=#I`Factorization eq 1 and I`Factorization[1][2] eq 1;
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
    // the direct test on S/I (index p^f, p in I): a product of two primes above p is reduced but
    // not a field (two fixed points of the Frobenius); a square containing p is not reduced
    for p in [2,3] do
        pp:=PrimesAbove(p*OneIdeal(E));
        if #pp ge 2 then
            assert not IsPrime(Ideal(E,ZBasis(pp[1]*pp[2])));
        end if;
        for P in pp do
            P2:=P*P;
            if Algebra(E)!p in P2 then
                assert not IsPrime(Ideal(E,ZBasis(P2)));
            end if;
        end for;
    end for;
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
    // IsMaximal without the maximal order: an order is maximal if and only if it has no singular
    // primes (computed from the F_p-algebras S/pS); the maximal order is not computed (at assertion
    // level 1: at level 2 SingularPrimes compares with the conductor)
    lvl:=GetAssertions();
    SetAssertions(1);
    for f in [ x^2+1, (x^2+1)*(x^2-3), x^4+16, (x^2+5)*(x^2+20), x^6+3*x^5+99*x^2+2997 ] do
        A:=EtaleAlgebra(f); E:=EquationOrder(A);
        b:=IsMaximal(E);
        // not computed if E is not maximal; E itself, stored as the maximal order, if it is
        assert (not b and not assigned A`MaximalOrder) or (b and IsIdentical(A`MaximalOrder,E));
        assert b eq (E eq MaximalOrder(A));
    end for;
    SetAssertions(lvl);
    printf ".";
    // a maximal order stored by IsMaximal works as the one built by MaximalOrder (its components
    // are declared maximal to Magma, see IsProductOfOrders): the same Picard and unit groups,
    // primes and factorizations
    for f in [ x^2+1, x^3-2 ] do
        A1:=EtaleAlgebra(f); E1:=EquationOrder(A1);
        assert IsMaximal(E1);
        A2:=EtaleAlgebra(f); O2:=MaximalOrder(A2);
        for O in [ E1, O2 ] do
            assert #PicardGroup(O) eq 1 and IsFinite(TorsionSubgroup(UnitGroup(O)));
        end for;
        assert { myHash(P) : P in PrimesAbove(6*OneIdeal(E1)) } eq { myHash(P) : P in PrimesAbove(6*OneIdeal(O2)) };
        assert { <myHash(e[1]),e[2]> : e in Factorization(30*OneIdeal(E1)) } eq { <myHash(e[1]),e[2]> : e in Factorization(30*OneIdeal(O2)) };
    end for;
    printf ".";
    // PrimesAboveRationalPrime: the same objects as PrimesAbove; for a maximal order, the primes
    // from the factorization in the components, with ramification index and inertia degree
    A:=EtaleAlgebra((x^2+5)*(x^3-2)); E:=EquationOrder(A); O:=MaximalOrder(A);
    for S in [ E, O ] do
        for p in [2,3,5] do
            pp:=PrimesAboveRationalPrime(S,p);
            assert forall{ P : P in PrimesAbove(p*OneIdeal(S)) | exists{ Q : Q in pp | IsIdentical(P,Q) } };
        end for;
    end for;
    assert forall{ P : P in PrimesAboveRationalPrime(O,2) | assigned P`RamificationIndex };
    // an order found maximal by IsMaximal: the primes of the maximal order built later from the
    // components are the objects found from S/pS; SingularPrimes also stores the maximal order
    lvl:=GetAssertions();
    SetAssertions(1);
    A:=EtaleAlgebra(x^3-2); E:=EquationOrder(A);
    assert IsMaximal(E) and IsIdentical(MaximalOrder(A),E);
    fp:=E`PrimesAboveRationalPrime[2];
    pl:=PlacesAboveRationalPrime(A,2); // built from the components: the objects found from E/2E
    assert forall{ P : P in pl | exists{ Q : Q in fp | IsIdentical(P,Q) } };
    A:=EtaleAlgebra(x^3-2); E:=EquationOrder(A);
    assert #SingularPrimes(E) eq 0 and assigned E`IsMaximal and E`IsMaximal and IsIdentical(A`MaximalOrder,E);
    SetAssertions(lvl);
    // IsBass stops at the first singular prime which is not Bass: the same answer as the definition
    A:=EtaleAlgebra((x^2+5)*(x^2+20)); E:=EquationOrder(A);
    for S in OverOrders(E) do
        T:=Order(ZBasis(S) : CheckIsKnownOrder:=false); // a copy, without stored singular primes
        assert IsBass(T) eq forall{ P : P in SingularPrimes(S) | IsBassAtPrime(S,P) };
    end for;
    // the factorizations and the primes above, passed by *, meet and + (see pass_primes in Idl.m),
    // are those of fresh copies
    A:=EtaleAlgebra(x^4+16); E:=EquationOrder(A);
    reg:=[ P : P in PrimesAbove(15*OneIdeal(E)) | IsInvertible(P) ]; // regular primes
    _:=[ Factorization(P) : P in reg ];
    X:=reg[1]*reg[2]; Y:=reg[2]*reg[3];
    fresh:=func< I | Ideal(E,ZBasis(I)) >;
    for Z in [ X*Y, X meet Y, X+Y ] do
        if Z ne OneIdeal(E) then
            assert assigned Z`Factorization;
            assert { <myHash(t[1]),t[2]> : t in Z`Factorization } eq { <myHash(t[1]),t[2]> : t in Factorization(fresh(Z)) };
        end if;
        assert { myHash(P) : P in PrimesAbove(Z) } eq { myHash(P) : P in PrimesAbove(fresh(Z)) };
    end for;
    sing:=SingularPrimes(E);
    for Z in [ sing[1]*reg[1], sing[1] meet reg[1], sing[1]+X ] do
        _:=PrimesAbove(sing[1]); _:=PrimesAbove(reg[1]); _:=PrimesAbove(X);
        assert { myHash(P) : P in PrimesAbove(Z) } eq { myHash(P) : P in PrimesAbove(fresh(Z)) };
    end for;
    // IsPrime and Factorization fill each other: a regular prime gets its factorization, and the
    // factorization decides IsPrime
    P:=fresh(reg[1]); P`IsInvertible:=true;
    assert IsPrime(P) and assigned P`Factorization and #P`Factorization eq 1;
    I:=fresh(X); _:=Factorization(I);
    assert assigned I`IsPrime and not I`IsPrime;
    printf ".";
    // p inert in S (S/pS a field), in fresh algebras: the prime is pS, generated by p
    for c in [ <x^2+1,3>, <x^2+1,7>, <x^3+x+1,2> ] do
        A:=EtaleAlgebra(c[1]); E:=EquationOrder(A);
        pp:=PrimesAbove(c[2]*OneIdeal(E));
        assert #pp eq 1 and pp[1] eq c[2]*OneIdeal(E) and IsPrime(pp[1]) and IsInvertible(pp[1]);
        assert not assigned A`MaximalOrder;
    end for;
    printf ".";
    // MinimalInteger of the primes: stored when they are built (from S/pS, from the components of
    // the maximal order, by PrimeBelow, by IsPrime), and read from the index for an ideal known to
    // be prime; the same as from the coordinates of 1
    A:=EtaleAlgebra((x^2+5)*(x^3-2)); E:=EquationOrder(A);
    fresh:=func< P | Ideal(Order(P),ZBasis(P)) >;
    pp:=&cat[ PrimesAbove(p*OneIdeal(E)) : p in [2,3,5] ]; // from E/pE
    O:=MaximalOrder(A);
    ppO:=&cat[ PrimesAbove(p*O) : p in [2,3,5] ]; // from the components
    ppB:=[ PrimeBelow(E,P) : P in ppO ];
    for P in pp cat ppO cat ppB do
        assert assigned P`MinimalInteger and P`MinimalInteger eq MinimalInteger(fresh(P));
        Q:=fresh(P); Q`IsPrime:=true; // marked prime: read from the index
        assert MinimalInteger(Q) eq MinimalInteger(fresh(P));
        Q:=fresh(P); assert IsPrime(Q) and assigned Q`MinimalInteger and Q`MinimalInteger eq P`MinimalInteger;
    end for;
    printf ".";
    SetAssertions(1);
    printf " all good!";

*/
