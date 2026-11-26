/////////////////////////////////////////////////////
// Copyright 2025.
// Stefano Marseglia, stefano.marseglia89@gmail.com
// https://stmar89.github.io/index.html
//
// Distributed under the terms of the CC-BY 4.0 licence.
// https://creativecommons.org/licenses/by/4.0/
/////////////////////////////////////////////////////



freeze;

declare verbose AlgEtQOrd, 3;

// the replacements of the working Z-basis of orders, ideals and modules, see
// replace_zbasis_if_smaller
declare verbose ZBasis, 2;

declare type AlgEtQOrd;

import "Elt.m" : CreateAlgEtQElt;

declare attributes AlgEtQ : MaximalOrder,
                           ProductOfEquationOrders,
                           EquationOrder,
                           KnownOrders; // An indexed set containing all the orders in the algebra
                                        // constructed so far.
                                        // We store them, because to avoid creating multiple copies
                                        // of the same order. In particular, this is more efficient
                                        // for the attributes. Looking up an order in it uses the
                                        // intrinsic Hash for orders. The attribute is populated at
                                        // creation of the order.

declare attributes AlgEtQOrd : IsMaximal,
                              IsProductOfOrders,
                              IsProductOfOrdersInFactorAlgebras,
                              Index,
                              ZBasis, // the working Z-basis, a record of format zbasis_fmt
                              CanonicalZBasis, // the LLL reduction of the HNF, a record of format
                                               // zbasis_fmt: see ZBasisLLLFromHNF
                              Generators,
                              OneIdeal,
                              Algebra,
                              TraceDualIdeal,
                              Conductor,
                              IsGorenstein,
                              IsBass,
                              Hash, // the sequence returned by myHash
                              HashValue; // the integer returned by Hash

///hide-all
intrinsic Print(A::AlgEtQOrd)
{Print the order.}
  printf"Order in %o", A`Algebra;
end intrinsic;

intrinsic IsCoercible(S::AlgEtQOrd, x::.) -> BoolElt, Any
{Return whether x is coercible into S and the result if so.}
    bool,x:=IsCoercible(Algebra(S),x);
    if bool and x in S then
        return true,x;
    else
        return false,_;
    end if;
end intrinsic;

//----------
// Auxiliary functions
//----------

crZQ:=function(M)
//given a matrix M over the Integers() changes the base ring to Rationals()
    return ChangeRing(M,Rationals());
end function;

crQZ:=function(M)
//given a matrix M over the Rationals() changes the base ring to Integers()
    return ChangeRing(M,Integers());
end function;

Columns:=function(M)
    return Rows(Transpose(M));
end function;

hnf:=function(M)
//input: a matrix over the Integers()
//output: the matrix in HNF, without the zero rows
// HermiteForm puts the zero rows at the bottom: we count them, instead of computing the rank (a
// second echelon form).
    H:=HermiteForm(M);
    r:=Nrows(H);
    while r gt 0 and IsZero(H[r]) do
        r-:=1;
    end while;
    if r lt Nrows(H) then
        H:=RowSubmatrix(H,r);
    end if;
    return H;
end function;

MatrixAtoQ:=function(gens)
//input:  a sequence of elements of an AlgEtQ A
//output: a matrix over Rationals() whose coordinates are the AbsoluteCoordinates of gens
    return Matrix([AbsoluteCoordinates(g) : g in gens]);
end function;

MatrixAtoZ:=function(gens)
//return a matrix with integers entries which rows the elements of gens
    return crQZ(MatrixAtoQ(gens));
end function;

MatrixQtoA:=function(A,P)
//given a matrix with Integer or Rational entries, returns a sequence of elements of A,
//corresponding to the rows.
// The absolute basis of A is the concatenation of the bases of the components of A (see
// AbsoluteBasis). So the j-th component of the element with absolute coordinates s is the element
// of the j-th number field whose coordinates are the j-th block of s: no arithmetic is needed.
    nf:=Components(A);
    assert Ncols(P) eq AbsoluteDimension(A);
    P:=ChangeRing(P,Rationals());
    degs:=[ Degree(K) : K in nf ];
    starts:=[ 1 + &+[ Integers() | degs[k] : k in [1..j-1] ] : j in [1..#nf] ];
    out:=[];
    for r in Rows(P) do
        s:=Eltseq(r);
        x:=CreateAlgEtQElt(A,< nf[j]!s[starts[j]..starts[j]+degs[j]-1] : j in [1..#nf] >);
        x`AbsoluteCoordinates:=s;
        Append(~out,x);
    end for;
    return out;
end function;

// The Hermite normal form of a full Z-lattice L in an étale algebra, and its hash.
// Let M be the matrix of the coordinates of a Z-basis of L, d a positive integer such that d*M is
// integral, and H the (upper triangular) Hermite normal form of d*M. The rational matrix (1/d)*H
// depends only on L: not on the Z-basis, nor on d, since the HNF of c*X is c times the HNF of X for
// c>0. Orders, fractional ideals and modules are hashed with hash_from_hnf(H,d), see myHash: so the
// hash can be stored for free wherever H is computed anyway, and conversely H can be read from the
// hash with hnf_from_hash. Both use only integer operations done by Magma's built-in functions.

lattice_hnf:=function(gens)
// Given elements of an étale algebra spanning a full Z-lattice L (possibly more than a Z-basis),
// returns H and d as above, where M is the matrix of the coordinates of gens and d its denominator.
    M:=MatrixAtoQ(gens);
    d:=Integers()!Denominator(M);
    return hnf(crQZ(d*M)),d;
end function;

hash_from_hnf:=function(H,d)
// Given H and d as above, returns the hash of L: the denominator e of (1/d)*H followed by the
// entries of the upper triangular part of (e/d)*H, row by row. With g the gcd of d and of the
// entries of H, e = d/g and (e/d)*H = H/g.
    dim:=Nrows(H);
    assert Ncols(H) eq dim; // a lattice of full rank: a degenerate one would get a truncated hash
    g:=GCD(d,GCD(Eltseq(H)));
    entries:=&cat[ Eltseq(H[i])[i..dim] : i in [1..dim] ];
    if g ne 1 then
        entries:=[ x div g : x in entries ];
    end if;
    return [d div g] cat entries;
end function;

hnf_from_hash:=function(hash,d)
// The inverse of hash_from_hnf: given the hash of L and an integer d>0 such that d*M is integral
// for the matrix M of a Z-basis of L, returns the HNF H of d*M. Since d*M is integral, d is a
// multiple of the denominator e=hash[1].
    e:=hash[1];
    assert d mod e eq 0;
    T:=UpperTriangularMatrix(hash[2..#hash]); // the entries are filled row by row, as in
                                              // hash_from_hnf
    return (d div e)*Matrix(Integers(),Nrows(T),Ncols(T),Eltseq(T));
end function;

// Z-bases of lattices: orders, fractional ideals and modules.
// A Z-basis of a lattice X is stored in a record of format zbasis_fmt: the working Z-basis in
// X`ZBasis (what ZBasis returns and ZBasisLLL reduces, as in v1.0.2), and the canonical one, the
// LLL reduction of the HNF of the lattice, in X`CanonicalZBasis (see ZBasisLLLFromHNF). The fields
// hold the same basis in different forms, and each one is computed only when it is needed:
// - A: the algebra where the elements live (the universe algebra for a module);
// - seq: the elements;
// - d, M: the integer matrix M = d*MatrixAtoQ(seq), where d is the denominator of MatrixAtoQ(seq).
//   It depends only on the lattice (it is the smallest d with d*X in ZZ^n): it is the same for all
//   the Z-bases of X, and it is X`Hash[1] when the hash is known;
// - Minv: the inclusion matrix MatrixAtoQ(seq)^-1 = d*M^-1, see inclusion_matrix;
// - lll: whether the basis is LLL-reduced;
// - hnf: whether M is the HNF of d*X.
// The records are written only by the functions below and by ZBasisLLL and ZBasisLLLFromHNF:
// the forms of a basis cannot get out of step (as the inclusion matrix of an old ZBasis could, when
// it was a separate attribute).
zbasis_fmt:=recformat< A, seq : SeqEnum, d : RngIntElt, M : Mtrx, Minv : Mtrx, lll : BoolElt,
                       hnf : BoolElt >;

zb_from_seq:=function(A,seq : lll:=false, hnf:=false)
// the record of the Z-basis seq of a lattice in A
    return rec< zbasis_fmt | A:=A, seq:=seq, lll:=lll, hnf:=hnf >;
end function;

zb_from_dM:=function(A,d,M : lll:=false, hnf:=false)
// the record of the Z-basis given by the rows of (1/d)*M, for an integer d>0 and an integer matrix
// M. Any such d is allowed: d and M are divided by their gcd, so that d is the denominator of
// (1/d)*M, as in zbasis_fmt.
    g:=GCD(d,GCD(Eltseq(M)));
    if g ne 1 then
        d:=d div g;
        M:=Matrix(Integers(),Nrows(M),Ncols(M),[ c div g : c in Eltseq(M) ]);
    end if;
    return rec< zbasis_fmt | A:=A, d:=d, M:=M, lll:=lll, hnf:=hnf >;
end function;

zb_seq:=function(X : att:="ZBasis")
// the elements of the Z-basis of X in the attribute att (which must be assigned), computed from d
// and M if needed
    Z:=X``att;
    if not assigned Z`seq then
        Z`seq:=MatrixQtoA(Z`A,(1/Z`d)*crZQ(Z`M));
        X``att:=Z;
    end if;
    return Z`seq;
end function;

zb_dM:=function(X : att:="ZBasis")
// d and M of the Z-basis of X in the attribute att (which must be assigned), see zbasis_fmt,
// computed from the elements if needed; d is read from the hash when it is known
    Z:=X``att;
    if not assigned Z`d then
        Q:=MatrixAtoQ(Z`seq);
        Z`d:=assigned X`Hash select X`Hash[1] else Integers()!Denominator(Q);
        Z`M:=crQZ(Z`d*Q);
        X``att:=Z;
    end if;
    return Z`d,Z`M;
end function;

zb_inv:=function(X : att:="ZBasis")
// the inclusion matrix of the Z-basis of X in the attribute att (which must be assigned)
    Z:=X``att;
    if not assigned Z`Minv then
        d,M:=zb_dM(X : att:=att);
        Z:=X``att; // zb_dM may have filled d and M
        Z`Minv:=d*crZQ(M)^-1;
        X``att:=Z;
    end if;
    return Z`Minv;
end function;

total_bits:=function(M)
// the number of bits of the entries of the integer matrix M: the measure used to compare two
// Z-bases of the same lattice (see replace_zbasis_if_smaller). A zero entry counts as one bit: the
// zeros of an HNF are cheap in the integer arithmetic.
    return &+[ Ilog2(Abs(c)+1)+1 : c in Eltseq(M) ];
end function;

replace_zbasis_if_smaller:=function(X,Z : Force:=false)
// Given a lattice X (an order, a fractional ideal or a module) and a record Z of format zbasis_fmt
// with another Z-basis of the same lattice, makes Z the working Z-basis of X if X has none, if
// Force is true, or if Z is smaller: its matrix has at most 3/4 of the bits of the matrix of the
// working Z-basis (see total_bits). The margin avoids replacing a basis, and with it the
// coordinates returned by AbsoluteCoordinates, for no measurable gain: in our measurements
// (2026-11-26) the operations on lattices take the same time with Z-bases whose sizes differ by
// less, and are slower only with much larger ones (e.g. a*Z-basis, or a skewed basis). Returns
// whether the working Z-basis was replaced. Every replacement of the working Z-basis goes through
// this function.
// The Generators of an ideal or a module with at least dim generators (e.g. the old Z-basis) are
// replaced with the new Z-basis, as ZBasisLLL does.
    if assigned X`ZBasis and not Force then
        d,M:=zb_dM(X);
        if not assigned Z`d then
            Z`d:=d;
            Z`M:=crQZ(d*MatrixAtoQ(Z`seq));
        end if;
        assert Z`d eq d; // the denominator depends only on the lattice
        bX:=total_bits(M);
        bZ:=total_bits(Z`M);
        if 4*bZ gt 3*bX then
            return false;
        end if;
        vprintf ZBasis, 1 : "replace_zbasis_if_smaller: the working Z-basis (%o bits) is replaced by a basis with %o bits\n", bX, bZ;
    end if;
    X`ZBasis:=Z;
    if Type(X) ne AlgEtQOrd and assigned X`Generators and #X`Generators ge AbsoluteDimension(Z`A) then
        X`Generators:=zb_seq(X);
    end if;
    return true;
end function;

sum_of_lattices_hnf:=function(Xs : extra:=[])
// Given orders or fractional ideals in the same algebra, returns H and d as in lattice_hnf for the
// lattice spanned by all of them (and by the elements in extra), from the integer matrices of their
// working Z-bases (see zb_dM): when they are known, the Z-bases are not converted to matrices
// again, and their elements are not needed.
    ds:=[]; Ms:=[* *]; // a list: the matrices of the records are square, the one of extra not
    if #extra gt 0 then
        Q:=MatrixAtoQ(extra);
        Append(~ds,Integers()!Denominator(Q));
        Append(~Ms,crQZ(ds[1]*Q));
    end if;
    for X in Xs do
        if not assigned X`ZBasis then
            _:=ZBasis(X);
        end if;
        dX,MX:=zb_dM(X);
        Append(~ds,dX);
        Append(~Ms,MX);
    end for;
    d:=LCM(ds);
    M:=(d div ds[1])*Ms[1];
    for k in [2..#Ms] do
        M:=VerticalJoin(M,(d div ds[k])*Ms[k]);
    end for;
    return hnf(M),d;
end function;

all_but_one:=function(seq, op, one)
// Given a sequence seq of n elements, an associative and commutative operation op (a function of
// two arguments, e.g. func< I,J | I meet J >) and its neutral element one, returns rests and suf,
// where rests[i] is op of all the seq[j] with j ne i, and suf[i] is op of seq[i..n] for i >= 2
// (suf[1] is not computed and is one, as suf[n+1]). They come from the prefixes and the suffixes:
// 3n-6 operations for n >= 2, instead of n(n-2) when each rests[i] is computed from scratch, and
// never more.
    n:=#seq;
    suf:=[ one : i in [1..n+1] ];
    for i in [n..2 by -1] do
        suf[i]:=i eq n select seq[n] else op(seq[i],suf[i+1]);
    end for;
    if n eq 1 then
        return [ one ],suf;
    end if;
    pre:=[ seq[1] ]; // pre[i] is op of seq[1..i], needed up to n-1
    for i in [2..n-1] do
        Append(~pre,op(pre[i-1],seq[i]));
    end for;
    rests:=[ i eq 1 select suf[2] else i eq n select pre[n-1] else op(pre[i-1],suf[i+1]) : i in [1..n] ];
    return rests,suf;
end function;

op_over_product:=function(lists, op)
// Given non-empty sequences lists[1..n] and a function op of two arguments, returns the sequence of
// op(...op(op(x_1,x_2),x_3)...,x_n) for all the tuples (x_1,...,x_n) in CartesianProduct(lists), in
// the same order (the last entry changes fastest). The result for a prefix (x_1,...,x_k) is
// computed once, and used by all the tuples starting with it: m_1*m_2 + m_1*m_2*m_3 + ... +
// m_1*...*m_n operations, for lists of sizes m_i, instead of (n-1)*m_1*...*m_n.
    n:=#lists;
    if n eq 1 then
        return lists[1];
    end if;
    out:=[];
    idx:=[ 1 : k in [1..n] ]; // the current tuple
    pre:=[ lists[1][1] ]; // pre[k] is the result for the first k entries, for k in [1..n-1]
    for k in [2..n-1] do
        Append(~pre,op(pre[k-1],lists[k][1]));
    end for;
    repeat
        Append(~out,op(pre[n-1],lists[n][idx[n]]));
        // the next tuple: the last index changes fastest
        k:=n;
        while k ge 1 and idx[k] eq #lists[k] do
            idx[k]:=1;
            k-:=1;
        end while;
        if k eq 0 then
            break;
        end if;
        idx[k]+:=1;
        // the results for the prefixes of length k, ..., n-1 change
        for j in [k..n-1] do
            pre[j]:=j eq 1 select lists[1][idx[1]] else op(pre[j-1],lists[j][idx[j]]);
        end for;
    until false;
    return out;
end function;

///hide-none

///# Orders in étale algebras over $\mathbb{Q}$
/// Let $A$ be an étale algebra over $\mathbb{Q}$. An `order` in $A$ is a subring $R$ of $A$ whose underlying additive group is a free abelian group of rank equal to the absolute dimension of $A$.
/// Given an order $R$ in $A$, we say that $S$ is an `overorder` of $R$ if $R\subseteq S$.
/// In MAGMA, orders have type `AlgEtQOrd`. Elements will always be considered as elements of the algebra.

///## Creation
/// Whenever we create an order, we populate the attributes `Generators` and `ZBasis`.
/// The `ZBasis` is put in a canonical form (we use the Hermite normal form of the numerators and divide by the denominators), unless the parameter `Check` is set to $0$ and the given elements are as many as the dimension of the algebra: then they are taken as the `ZBasis`. The entries of the corresponding upper triangular matrix uniquely determine the order and are used for hashing it.

create_order:=function(A,gens,zb,hash,CheckIsKnownOrder)
// Creates the order of A with Z-basis zb, generated by gens (when gens has more than dim elements,
// the generators are the Z-basis, see Generators), with the given hash if hash is not empty (see
// hash_from_hnf). If CheckIsKnownOrder is true and the order is already known, returns the known
// one. The stored Z-basis is LLL-reduced.
// If zb is empty, the hash must be given: the Z-basis is then the LLL reduction of the HNF read
// from the hash, computed only if the order is not already known, and its elements only when they
// are asked for (see zbasis_fmt). This avoids converting the HNF basis to elements and back
// (MatrixQtoA, MatrixAtoQ and the denominator of the matrix cost about 6 times the LLL itself in
// OverOrders_3.12).
    dim:=AbsoluteDimension(A);
    R:=New(AlgEtQOrd);
    R`Algebra:=A;
    if #hash gt 0 then
        R`Hash:=hash;
    end if;
    if #zb eq 0 then
        assert #hash gt 0;
        R`ZBasis:=zb_from_dM(A,hash[1],hnf_from_hash(hash,hash[1]) : hnf:=true);
    else
        assert #zb eq dim;
        R`ZBasis:=zb_from_seq(A,zb);
    end if;
    if CheckIsKnownOrder then
        R0:=R;
        IsKnownOrder(~R);
        if not IsIdentical(R,R0) then // R is the stored copy of a known order
            return R;
        end if;
    end if;
    if #gens gt 0 and #gens le dim then // otherwise Generators returns the Z-basis
        R`Generators:=gens;
    end if;
    assert2 One(A) in R;
    ZBasisLLL(R); // for better coefficients: it uses the HNF stored in the hash, if known
    return R;
end function;

order_from_sum:=function(Xs)
// The order spanned as a Z-module by the orders and fractional ideals in Xs, whose sum is assumed
// to be closed under multiplication, as Order(&cat[ ZBasis(X) : X in Xs ] : Check:=0) (the same
// output), but from the integer matrices of the Z-bases (see sum_of_lattices_hnf). Examples: S+T
// for overorders S and T of R which differ from R at distinct primes, S+P^kO.
    A:=Algebra(Xs[1]);
    H,d:=sum_of_lattices_hnf(Xs);
    assert Nrows(H) eq AbsoluteDimension(A);
    hash:=hash_from_hnf(H,d);
    if GetAssertions() ge 2 then // the lattice is closed under multiplication
        zbH:=MatrixQtoA(A,(1/d)*crZQ(H));
        H2,d2:=lattice_hnf(zbH cat [ z1*z2 : z1,z2 in zbH ]);
        assert hash_from_hnf(H2,d2) eq hash;
    end if;
    return create_order(A,[],[],hash,true);
end function;

/// Returns the order generated as an algebra over $\mathbb{Z}$ by a given sequence of elements. The parameter `Check` (default $100$) determines how many times the program tries to obtain a multiplicatively closed lattice by adding the product of the given generators. If `Check` is $0$ then this step is skipped: the $\mathbb{Z}$-span of the given elements is assumed to be closed under multiplication, which is not verified, and it is returned as an order. The parameter `CheckIsKnownOrder` determines whether the program checks if the order is already known, i.e. in the attribute `KnownOrders` of the algebra. The default value is true. The stored $\mathbb{Z}$-basis is LLL reduced.
intrinsic Order( gens::SeqEnum[AlgEtQElt] : Check:=100 , CheckIsKnownOrder:=true ) -> AlgEtQOrd
{Returns the order generated as an algebra over ZZ by a given sequence of elements. The parameter Check (default 100) determines how many times the program tries to obtain a multiplicatively closed lattice by adding the product of the given generators. If Check is 0 then this step is skipped: the ZZ-span of the given elements is assumed to be closed under multiplication, which is not verified, and it is returned as an order. The vararg CheckIsKnownOrder determines whether the program checks if the order is already known, i.e. in the attribute KnownOrders of the algebra. The default value is true. The stored ZBasis is LLL-reduced.}
    A:=Universe(gens);
    dim:=AbsoluteDimension(A);
    // We remove duplicates, which matter when the generators are used as they are: in the closure
    // loop, and when they are taken as a Z-basis. With Check eq 0 and more than dim generators only
    // their Z-span (an HNF) is used, so the duplicates do not matter, and removing them costs a
    // comparison of elements per pair with the same hash (about 25% of OverOrders for many
    // overorders).
    if Check gt 0 or #gens le dim then
        gens:=Setseq(Seqset(gens));
    end if;
    hash:=[]; // the hash of the order, when its HNF is computed below
    if Check gt 0 then
        require forall{g : g in gens | IsIntegral(g) } : "The elements are not all integral, hence do not generate an order.";
    else
        // the Z-span of the elements is assumed to be an order: in particular, they are integral
        assert2 forall{g : g in gens | IsIntegral(g) };
    end if;
    if Check gt 0 then
        gens0:=gens;
        if not One(A) in gens then
            gens0:=[One(A)] cat gens0;
        end if;
        M:=MatrixAtoQ(gens0);
        d:=Denominator(M);
        P:=crZQ(hnf(crQZ(d*M)))/d;
        gens0:=MatrixQtoA(A,P);
        fail:=0;
        go:=false;
        repeat
            P_old:=P;
            // the algebra is commutative: the products g_i*g_j with i le j suffice
            gens0 cat:= &cat[ [ gens0[i]*gens0[j] : j in [i..#gens0] ] : i in [1..#gens0] ];
            M:=MatrixAtoQ(gens0);
            d:=Denominator(M);
            M:=crQZ(M*d);
            P:=crZQ(hnf(M))/d;
            gens0:=MatrixQtoA(A,P);
            fail +:=1;
            if P cmpeq P_old then
                if Rank(P_old) eq dim then
                    go:=true;
                else
                    error "The algebra over ZZ generated by the given set is not a full ZZ-lattice in the étale algebra.";
                end if;
            end if;
        until go or fail gt Check;
        if not go then
            error "The program failed to compute the Z-algebra generated by the given set. Try to increase the parameter Check.";
        end if;
        e:=Integers()!Denominator(P);
        hash:=hash_from_hnf(crQZ(e*P),e);
        zb:=[]; // the LLL-reduced Z-basis is computed from the hash in create_order
    elif #gens gt dim then
    // We assume that gens generates over Z a multiplicatively closed lattice: we only need its
    // Z-basis (one HNF). Examples:
    // - S_1+...+S_n, for overorders S_i of an order S which differ from S at pairwise distinct
    //   primes of S: at each prime of S the sum is locally equal to one of the S_i, hence it is
    //   closed under multiplication;
    // - S+I, for an order S and an integral ideal I of an overorder of S.
        H,d:=lattice_hnf(gens);
        hash:=hash_from_hnf(H,d);
        zb:=[]; // the LLL-reduced Z-basis is computed from the hash in create_order
        // the lattice is closed under multiplication: it has the same hash as the lattice spanned
        // with the products
        if GetAssertions() ge 2 then
            zbH:=MatrixQtoA(A,(1/d)*crZQ(H));
            H2,d2:=lattice_hnf(zbH cat [ z1*z2 : z1,z2 in zbH ]);
            assert hash_from_hnf(H2,d2) eq hash;
        end if;
    else
    //we assume that gens is a ZBasis of a multiplicativly closed lattice
        require #gens eq AbsoluteDimension(A) : "with Check:=0 at least as many elements as the dimension are needed";
        zb:=gens;
        // as in the branch above: a lattice of full rank, closed under multiplication
        if GetAssertions() ge 2 then
            H,d:=lattice_hnf(zb);
            H2,d2:=lattice_hnf(zb cat [ z1*z2 : z1,z2 in zb ]);
            assert hash_from_hnf(H2,d2) eq hash_from_hnf(H,d);
        end if;
    end if;

    return create_order(A,gens,zb,hash,CheckIsKnownOrder);
end intrinsic;

intrinsic Order(A::AlgEtQ , orders::Tup) -> AlgEtQOrd
{Given an étale algebra and a tuple of orders in the components of the étale algebra, it generates the direct sum order.}
    nf,embs:=Components(A);
    // ISA, not Type(S) in {...}: the orders of cyclotomic and quadratic fields have the types
    // RngCyc and RngQuad, which are subtypes of RngOrd.
    require #nf eq #orders and forall{ i : i in [1..#orders] | NumberField(orders[i]) eq nf[i]} and
forall{ S : S in orders | ISA(Type(S),RngOrd) } : "The second input should be a sequence of orders in the number fields defining A";
    gens:=&cat[ [ embs[i](z) : z in AbsoluteBasis(orders[i]) ] : i in [1..#nf] ];
    O:=Order(gens : Check:=0); // the gens given generate a multiplicatively closed lattice in A
    assert2 O eq Order(gens);
    O`IsProductOfOrders:=<true,orders>;
    return O;
end intrinsic;

///## Equation orders
/// If $a$ is the stored primitive element of $A$ then we refer to the order $\mathbb{Z}[a]$ as the `equation order` of $A$.

intrinsic EquationOrder(A::AlgEtQ) -> AlgEtQOrd
{Returns the monogenic order of the étale algebra, which depends on the stored primitive element.}
    if not assigned A`EquationOrder then
        require PrimeField(A) eq BaseField(A) : "Defined only for algebras over the Rationals()";
        require IsIntegral(PrimitiveElement(A)) : "the primitive element is not integral";
        pow:=PowerBasis(A);
        E:=Order(pow : Check:=0); //the PowerBasis generates a multiplicatively closed lattice
        assert2 E eq Order([pow[2]]);
        E`Generators:=[pow[2]]; // I want the generator set to be small, so I put only the primitive
                                // element
        A`EquationOrder:=E;
    end if;
    return A`EquationOrder;
end intrinsic;

intrinsic ProductOfEquationOrders(A::AlgEtQ)->AlgEtQOrd
{Given an étale algebra A, returns the order consisting of the product of the equation orders of the number fields.}
    if not assigned A`ProductOfEquationOrders then
        E:=Order( A , <EquationOrder(E) : E in Components(A)> );
        A`ProductOfEquationOrders := E;
    end if;
    return A`ProductOfEquationOrders ;
end intrinsic;

///## Maximal Order
/// If $A$ has components $K_1,\ldots,K_n$ then product $\mathcal{O}_A$ of the ring of integers $\mathcal{O}_{K_i}$ of $K_i$ is an order, which contains every other order of $A$. For this reason, we refer to $\mathcal{O}_A$ as the `maximal order` of $A$.

intrinsic MaximalOrder(A::AlgEtQ)->AlgEtQOrd
{Returns the maximal order of the étale algebra A. It is the direct sum of the ring of integers of the number fields composing the algebra.}
    if not assigned A`MaximalOrder then
        O:=Order( A , <MaximalOrder(E) : E in Components(A)> );
        O`IsMaximal:=true;
        A`MaximalOrder:=O;
    end if;
    return A`MaximalOrder;
end intrinsic;

intrinsic IsMaximal(S::AlgEtQOrd) -> BoolElt
{Returns whether the given order is the maximal order of the étale algebra.}
    if not assigned S`IsMaximal then
        O:=MaximalOrder(Algebra(S));
        S`IsMaximal:=S eq O;
    end if;
    return S`IsMaximal;
end intrinsic;

///## The attribute `KnownOrders`
/// An étale algebra $A$ has an attribute `KnownOrders`, in which all previously constructed orders are stored.
/// Whenever the creation intrinsic for orders is called, the program checks if the orders has already been constructed. This is done to prevent the creation of the same order multiple times. This check can be turned off with the appropriate parameter.

intrinsic IsKnownOrder(~R::AlgEtQOrd)
{This procedure checks whether the given order is already in the list of known orders of the parent algebra. If so then it replaces the given order with the stored copy. If not, it adds is added to the storage.}
    A:=Algebra(R);
    // A`KnownOrders is an indexed set: looking up R uses the intrinsic Hash for orders.
    if not assigned A`KnownOrders then
        A`KnownOrders:={@ R @};
    else
        i:=Index(A`KnownOrders,R);
        if i gt 0 then
            R:=A`KnownOrders[i];
        else
            // We take the indexed set out of the attribute, so that Include modifies it in place:
            // otherwise the whole set is copied, and every order in it hashed again, at every new
            // order.
            K:=A`KnownOrders;
            delete A`KnownOrders;
            Include(~K,R);
            A`KnownOrders:=K;
        end if;
    end if;
end intrinsic;

///## Main attributes of an order

intrinsic Algebra(S::AlgEtQOrd) -> AlgEtQ
{Returns the algebra of the order.}
  return S`Algebra;
end intrinsic;

/// Returns the stored $\mathbb{Z}$-basis of the order, the working one: it is replaced by `ZBasisLLL`, and also when a smaller $\mathbb{Z}$-basis of the same lattice becomes known (for example its Hermite normal form, computed by `myHash`), see `ReplaceZBasis`. So the coordinates returned by `AbsoluteCoordinates` refer to the $\mathbb{Z}$-basis stored at the time of the call.
intrinsic ZBasis(S::AlgEtQOrd)->SeqEnum[AlgEtQElt]
{Return a Z-basis of the order, the stored working one: it is replaced by ZBasisLLL, and also when a smaller Z-basis of the same lattice becomes known (for example its Hermite normal form, computed by myHash), see ReplaceZBasis. So the coordinates returned by AbsoluteCoordinates refer to the Z-basis stored at the time of the call.}
    // the record is assigned at creation, its elements possibly only now (see zbasis_fmt)
    return zb_seq(S);
end intrinsic;

/// Return a set of generators as a $\mathbb{Z}$-algebra of the order.
intrinsic Generators(S::AlgEtQOrd)->SeqEnum[AlgEtQElt]
{Return a set of generators (as a Z-algebra) of the order.}
    // assigned at creation, unless the order was created from more than dim elements: then the
    // generators are the Z-basis
    if not assigned S`Generators then
        S`Generators:=ZBasis(S);
    end if;
    return S`Generators;
end intrinsic;

///## Equality
/// Equality of orders is performed using the `Hash` attribute which is constructed as follows.
/// Let $S$ be an order.
/// Let $P$ be the upper triangular Hermite normal form of the integer square matrix $d\cdot M$ where $M$ is the matrix whose rows are the coefficients of a $\mathbb{Z}$-basis of $S$ and $d$ is the least common denominator of its entries.
/// The `Hash` of $S$ is defined to be the sequence consisting of the least common denominator of $\frac{1}{d}\cdot P$ and the entries of the upper triangular part of $\frac{1}{d}\cdot P$.
/// This hashing method has no collisions and it is independent of the choice of $\mathbb{Z}$-basis from which we start the procedure.
/// We observed that applying the in-built Hash function to the sequence defined above, while giving a smaller hash, it often lead to collisions. For this reason, equality is tested using the whole sequence.
/// The intrinsic `Hash` returns an integer computed from this sequence, which is used by Magma to build sets, indexed sets and associative arrays of orders: only orders with the same hash value are compared with `eq`, so collisions only cost an extra comparison.

intrinsic myHash(S::AlgEtQOrd)->SeqEnum[RngIntElt]
{Hash function for AlgEtQOrd.}
// two orders have the hash iff they are equal
    if not assigned S`Hash then
        d,M:=zb_dM(S);
        H:=hnf(M);
        assert2 IsUpperTriangular(H);
        S`Hash:=hash_from_hnf(H,d);
        // the HNF is a Z-basis of S: it becomes the working one if it is smaller
        _:=replace_zbasis_if_smaller(S,zb_from_dM(Algebra(S),d,H : hnf:=true));
    end if;
    return S`Hash;
end intrinsic;

/// Returns a hash value of the order, so that equal orders have the same hash value. It is computed from the sequence returned by `myHash` and it is stored in the attribute `HashValue`. Magma uses it to build sets, indexed sets and associative arrays of orders: only orders with the same hash value are compared with `eq`. Different orders might have the same hash value: this only costs an extra comparison.
intrinsic Hash(S::AlgEtQOrd) -> RngIntElt
{Returns a hash value of the order, so that equal orders have the same hash value. It is computed from the sequence returned by myHash and it is stored in the attribute HashValue. Magma uses it to build sets, indexed sets and associative arrays of orders: only orders with the same hash value are compared with eq. Different orders might have the same hash value: this only costs an extra comparison.}
    if not assigned S`HashValue then
        // The universe matters, e.g. Hash([1,2]) ne Hash([1/1,2/1]), so we make sure it is
        // Integers().
        S`HashValue:=Hash(ChangeUniverse(myHash(S),Integers()));
    end if;
    return S`HashValue;
end intrinsic;

intrinsic 'eq'(O1::AlgEtQOrd,O2::AlgEtQOrd)->BoolElt
{Checks equality of orders.}
    if IsIdentical(O1,O2) then
        return true;
    end if;
    require Algebra(O1) cmpeq Algebra(O2) : "the orders must be defined in the same etale algebra";
    if
    //O1`Generators eq O2`Generators or O1`ZBasis eq O2`ZBasis or //this seems to be slowe
    myHash(O1) eq myHash(O2) then
        out:=true;
    else
        out:=false;
    end if;
    vprintf AlgEtQOrd,3 : "output of equality = %o\n",out;
    return out;
end intrinsic;

///## Inclusion and coordinates
/// Inclusion testing of elements, orders and ideals in a fixed order is performed by multiplying by the `inclusion matrix`. This matrix, which is stored in an attribute, is the inverse of the matrix with coefficients the $\mathbb{Z}$-basis of the order. If the output of the multiplication has integer coefficients then we have an inclusion.
/// The same matrix can be used to obtain the coordinates of a sequence of elements with respect to the $\mathbb{Z}$-basis.

inclusion_matrix:=function(O)
// The inclusion matrix of the working Z-basis of O (an order, a fractional ideal or a module): the
// inverse of the matrix MatrixAtoQ(ZBasis(O)), stored in the record of the Z-basis (see
// zbasis_fmt), so that it always belongs to the current Z-basis.
// This allows to check if an element is in O by simply making one matrix multiplication, rather
// then solving a system (as in AbsoluteCoordinates). It is much faster.
    if not assigned O`ZBasis then
        _:=ZBasis(O); // an ideal or a module whose Z-basis is not computed yet
    end if;
    return zb_inv(O);
end function;

all_in:=function(seq,X)
// Given a sequence of elements of the algebra of X, an order, a fractional ideal or a module
// (anything with a ZBasis), returns whether all of them are in X: one matrix for all the elements,
// one product with the inclusion matrix, and one test that the product is integral (CanChangeRing
// is about 8 times faster than testing the entries one by one with IsCoercible).
    if #seq eq 0 then
        return true;
    end if;
    Minv:=inclusion_matrix(X);
    is_in:=CanChangeRing(MatrixAtoQ(seq)*Minv,Integers()); // only the first of its two values
    return is_in;
end function;

intrinsic 'in'(x::AlgEtQElt,O::AlgEtQOrd) -> BoolElt
{Inclusion of an element in the order.}
    require Algebra(x) cmpeq Algebra(O) : "the algebra is not the same";
    /* //this check succeeds very rarely, and makes the code a bit slower.
    if (assigned O`Generators and x in O`Generators)
        or (assigned O`ZBasis and x in O`ZBasis) then
        return true;
    end if;
    */
    return all_in([x],O);
end intrinsic;

intrinsic 'in'(x::RngIntElt,O::AlgEtQOrd) -> BoolElt
{Inclusion of elements.}
    return (Algebra(O)!x) in O;
end intrinsic;

intrinsic 'in'(x::FldRatElt,O::AlgEtQOrd) -> BoolElt
{Inclusion of elements.}
    return (Algebra(O)!x) in O;
end intrinsic;

intrinsic 'subset'(O1 :: AlgEtQOrd, O2 :: AlgEtQOrd) -> BoolElt
{Checks if the first argument is inside the second.}
  require Algebra(O1) cmpeq Algebra(O2) : "The orders must be in the same algebra.";
  if not Index(O2, O1) in Integers() then
    return false;
  end if;
  return all_in(Generators(O1),O2); // ring generators suffice: O2 is closed under products
end intrinsic;

/// Given an element $x$ and an order $S$, returns the coordinates of the element with respect to the stored $\mathbb{Z}$-basis of $S$.
intrinsic AbsoluteCoordinates(x::AlgEtQElt,S::AlgEtQOrd) -> SeqEnum
{Given an element x and an order S, returns the coordinates of the element with respect to the stored Z-basis of S.}
    return AbsoluteCoordinates([x],S)[1];
end intrinsic;

/// Returns the absolute coordinates of a sequence of elements  with respect to the stored $\mathbb{Z}$-basis of the given order.
intrinsic AbsoluteCoordinates(seq::SeqEnum[AlgEtQElt],O::AlgEtQOrd) -> SeqEnum
{Returns the absolute coordinates of a sequence of elements  with respect to the stored Z-basis of the given order.}
    require forall{ x : x in seq | Algebra(x) cmpeq Algebra(O)} : "the algebra is not the same";
    Minv:=inclusion_matrix(O);
    M:=MatrixAtoQ(seq)*Minv;
    return [Eltseq(r) : r in Rows(M)];
end intrinsic;

///## Special and random elements in an order

intrinsic One(S::AlgEtQOrd)->AlgEtQElt
{Returns the multiplicative unit of the order.}
    A:=Algebra(S);
    return One(A);
end intrinsic;

intrinsic Zero(S::AlgEtQOrd)->AlgEtQElt
{Returns the additive neutral element of the order.}
    A:=Algebra(S);
    return Zero(A);
end intrinsic;

/// Returns a random element of the given order with coefficients with respect to the stored $\mathbb{Z}$-basis bounded by `bd`. One can allow zero-divisors using the parameter `ZeroDivisorsAllowed`, which is set to false by default.
intrinsic Random(O::AlgEtQOrd , bd::RngIntElt : ZeroDivisorsAllowed:=false ) -> AlgEtQElt
{Returns a random element of O with coefficcients with respect to the stored Z-basis bounded by bd. One can allow zero-divisors using the parameter ZeroDivisorsAllowed, which is set to false by default.}
    require bd gt 0 : "The bound needs to be a positive integer.";
    B := ZBasis(O);
    if ZeroDivisorsAllowed then
       elt:=DotProduct([ Random([-bd..bd]) : b in B ],B);
    else
        repeat
            elt:=DotProduct([ Random([-bd..bd]) : b in B ],B);
        until not IsZeroDivisor(elt);
    end if;
    return elt;
end intrinsic;

/// Returns a random element of the given order. The coefficients with respect to the stored $\mathbb{Z}$-basis can be bounded by setting the parameter `CoeffRange`. One can allow zero-divisors using the parameter `ZeroDivisorsAllowed`, which is set to false by default.
intrinsic Random(O::AlgEtQOrd : CoeffRange:=3, ZeroDivisorsAllowed:=false ) -> AlgEtQElt
{Returns a random (small coefficient) element of O. The range of the random coefficients can be increased by giving the optional argument CoeffRange. One can allow zero-divisors using the optional argument "ZeroDivisorsAllowed", which is set to false by default.}
      return Random(O,CoeffRange : ZeroDivisorsAllowed:=ZeroDivisorsAllowed);
end intrinsic;


//----------
// Index
//----------

/// Given an order $T$ returns the determinant of the change of basis from a $\mathbb{Z}$-basis of $T$ to the basis of the parent algebra.
intrinsic Index(T::AlgEtQOrd) -> FldRatElt
{Given an order T returns the determinant of the change of basis from a Z-basis of T to the basis of the parent algebra.}
  if not assigned T`Index then
    // the determinant of (1/d)*M, see zbasis_fmt
    d,M:=zb_dM(T);
    T`Index := Abs(Determinant(M))/d^Nrows(M);
  end if;
  return T`Index;
end intrinsic;

/// Given two orders $S$ and $T$ returns the ratio of their indices. If $T subset S$ the outputs is the index of the inclusion, that is, it equals $[S:T] = \#S/T$.
intrinsic Index(S::AlgEtQOrd, T::AlgEtQOrd) -> FldRatElt
{Given two orders S and T returns the ratio of their indices. If T subset S the outputs is the index of the inclusion, that is, it equals [S:T] = \#S/T.}
  require Algebra(T) cmpeq Algebra(S) : "the orders must be in the same algebra";
  elt := Index(T)/Index(S);
  if IsCoercible(Integers(), elt) then
    elt := Integers() ! elt;
  end if;
  return elt;
end intrinsic;

///## Product and intersection of orders

intrinsic '*'(O1::AlgEtQOrd,O2::AlgEtQOrd)->AlgEtQOrd
{Returns the order generated by the orders O1 and O2.}
    require Algebra(O1) cmpeq Algebra(O2) : "the orders must be defined in the same algebra";
    if O1 eq O2 then
        return O1;
    end if;
    if assigned O1`IsMaximal and O1`IsMaximal then
    //calling IsMaximal(O1) would trigger the computation of the maximal order, which might be
    //expensive
        return O1;
    end if;
    if assigned O2`IsMaximal and O2`IsMaximal then
        return O2;
    end if;
    if O1 subset O2 then
        return O2;
    end if;
    if O2 subset O1 then
        return O1;
    end if;
    gens:=Setseq( Seqset(Generators(O1) cat Generators(O2)) );
    return Order(gens);
end intrinsic;

meet_zbasis:=function(X1,X2)
// Given two orders or fractional ideals (or one of each) in the same algebra, returns the record
// (see zbasis_fmt) of the HNF Z-basis of the intersection of their lattices, without its elements,
// and its hash.
// The intersection depends only on the lattices: when the hash of X is known, d times the HNF of
// its lattice is read from the hash (hnf_from_hash, for any multiple d of hash[1]); otherwise the
// integer matrix of its working Z-basis is used (see zb_dM).
    A:=Algebra(X1);
    N:=AbsoluteDimension(A);
    Ms:=[* *]; ds:=[];
    for X in [* X1, X2 *] do
        if assigned X`Hash then
            Append(~Ms,X`Hash);
            Append(~ds,X`Hash[1]);
        else
            if not assigned X`ZBasis then
                _:=ZBasis(X);
            end if;
            dX,MX:=zb_dM(X);
            Append(~Ms,MX);
            Append(~ds,dX);
        end if;
    end for;
    d:=LCM(ds);
    scaled:=func< M,dX | Type(M) eq SeqEnum select hnf_from_hash(M,d) else (d div dX)*M >;
    MI:=scaled(Ms[1],ds[1]);
    MJ:=scaled(Ms[2],ds[2]);
    MZ:=ZeroMatrix(Integers(),N,N);
    M:=Matrix(HorizontalJoin(VerticalJoin(MI,MJ),VerticalJoin(MI,MZ)));
    M:=hnf(M);
    //the hnf of zb1 meet zb2 is the lower-right quadrant of M
    P:=Matrix(N,N,[M[i,j] : i,j in [N+1..2*N]]);
    // A block of a matrix in HNF is in HNF: P is the HNF of d times the matrix of a Z-basis of the
    // intersection, so its hash is free (see lattice_hnf).
    return zb_from_dM(A,d,P : hnf:=true),hash_from_hnf(P,d);
end function;

intrinsic 'meet'(O1::AlgEtQOrd,O2::AlgEtQOrd)->AlgEtQOrd
{Returns the order given by the intersection of two given orders.}
    require Algebra(O1) cmpeq Algebra(O2) : "the orders must be defined in the same algebra";
    if assigned O1`IsMaximal and O1`IsMaximal then
    //calling IsMaximal(O1) would trigger the computation of the maximal order, which might be
    //expensive
        return O2;
    end if;
    if assigned O2`IsMaximal and O2`IsMaximal then
        return O1;
    end if;
    if O1 eq O2 then
        return O1;
    end if;
    if O1 subset O2 then
        return O1;
    end if;
    if O2 subset O1 then
        return O2;
    end if;
    _,hash:=meet_zbasis(O1,O2);
    // The intersection of two orders is multiplicatively closed. Its hash is known: the order is
    // created from it (see create_order), and it is looked up among the known orders.
    O:=create_order(Algebra(O1),[],[],hash,true);
    assert2 O eq Order(ZBasis(O));
    return O;
end intrinsic;

///## Orders which are direct products of orders

intrinsic IsProductOfOrdersInComponents(O::AlgEtQOrd)->BoolElt, Tup
{Returns whether the argument is a product of orders in the components of its parent algebra. If so, it returns also a tuple containing these orders.}
    return IsProductOfOrders(O);
end intrinsic;

///ditto
intrinsic IsProductOfOrders(O::AlgEtQOrd)->BoolElt, Tup
{ " } // "
    if not assigned O`IsProductOfOrders then
        A:=Algebra(O);
        idem:=OrthogonalIdempotents(A);
        test:=all_in(idem,O);
        if test then
            nf,_,projs:=Components(A);
            if #nf eq 1 then
                orders:=< Order([ projs[1](z) : z in ZBasis(O) ])  >;
            else
                zb:=ZBasis(O);
                require PrimeField(A) eq BaseField(A) : "We require the algebra to be defined over the Rationals(). The reason is that we don't have a functioning method to create orders in relative number fields which are not over the maximal order of the base field. One can still check if the order O is a product by checking if OrthogonalIdempotents(A) are in O.";
                orders:=< Order( [ projs[i](z) : z in zb ] ) : i in [1..#idem] >;
                assert2 O eq Order(A,orders);
            end if;
            O`IsProductOfOrders:=<true,orders>;
        else
            O`IsProductOfOrders:=<false,<>>;
        end if;
    end if;
    return Explode(O`IsProductOfOrders);
end intrinsic;

intrinsic IsProductOfOrdersInFactorAlgebras(S::AlgEtQOrd)->BoolElt,SeqEnum[AlgEtQElt]
{Returns whether the given order is a product of orders living in some factor algebras of the parent algebra.
This is equivalent for the given order to contain some idempotents of the algebra other than 0 and 1. If this is the case, it returns also the idempotents.}
    if not assigned S`IsProductOfOrdersInFactorAlgebras then
        A:=Algebra(S);
        idem:=[ i : i in Idempotents(A) | i in S]; // idempotents in S
        output:=[PowerStructure(AlgEtQElt)|];
        if Seqset(idem) eq { One(A), Zero(A) } then
            S`IsProductOfOrdersInFactorAlgebras:=<false,output>;
        else
            test,D:=IsPowerOf(#idem,2);
            assert test;
            F2:=GF(2);
            V:=VectorSpace(F2,#Components(A));
            idem_inV:=[ V![ c eq 1 select F2!1 else F2!0 : c in Components(i) ] : i in idem ];
            arr:=AssociativeArray();
            for j in [1..#idem] do
                i:=idem[j];
                v:=idem_inV[j];
                n:=#[ x : x in Eltseq(v) | x eq 1 ];
                if not IsDefined(arr,n) then
                    arr[n]:=[];
                end if;
                Append(~arr[n],j);
            end for;
            W:=sub<V|>;
            n:=1;
            while Dimension(W) lt D do
                if IsDefined(arr,n) then
                    for j in arr[n] do
                        v:=idem_inV[j];
                        if not v in W then
                            Append(~output,idem[j]);
                            W:=sub<V| W,v>;
                            if Dimension(W) eq D then
                                break;
                            end if;
                        end if;
                    end for;
                end if;
                n+:=1;
            end while;
            assert #output eq D;
            assert2 OneIdeal(S) eq Ideal(S , [ z*g : z in ZBasis(S) , g in output ]);

            S`IsProductOfOrdersInFactorAlgebras:=<true, output>;
        end if;
    end if;
    return Explode(S`IsProductOfOrdersInFactorAlgebras);
end intrinsic;

///# Example 3
/// ```
/// _<x>:=PolynomialRing(Integers());
/// // We consider the following three number fields
/// K1:=NumberField(x^2-2);
/// K2:=NumberField(x^2-3);
/// K3:=NumberField(x^2-5);
/// // We define the product étale algebra A and the factor algebras B and C consisting of only the first two components and the last one, respectively.
/// B:=EtaleAlgebra([K1,K2]);
/// C:=EtaleAlgebra([K3]);
/// A,embs,projs:=DirectProduct([B,C]);
/// // The maximal order of A is the product of the three ring of integers.
/// OA:=MaximalOrder(A);
/// IsProductOfOrdersInComponents(OA);
/// // The equation order of A is not a product in any factor algebra
/// EA:=EquationOrder(A);
/// IsProductOfOrdersInFactorAlgebras(EA);
/// // Now we construct an order that is a product of an order in B and one in C, but does not admit further splittings.
/// a:=PrimitiveElement(A);
/// e1:=A![1,0,0];
/// e2:=A![0,1,1];
/// R:=Order([a*e1,a*e2]);
/// IsProductOfOrdersInFactorAlgebras(R);
/// ```

/* TESTS

    printf "### Testing Orders:";
    //AttachSpec("../spec");
    SetVerbose("AlgEtQOrd",1);
    SetAssertions(2);

    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    O1:=Order(Basis(A));
    _:=ZBasis(O1);
    _:=Generators(O1);
    O2:=Order(AbsoluteBasis(A) : Check:=0);
    O3:=Order(AbsoluteBasis(A));
    _:=O1 eq O2;
    _:=O2 eq O3;
    assert EquationOrder(A) ne ProductOfEquationOrders(A);
    printf ".";

    // we test KnownOrders
    OA:=MaximalOrder(A);
    O:=Order(ZBasis(OA));
    assert assigned O`IsMaximal; //at creation O<-OA since we call IsKnownOrder
    assert O eq OA;
    printf ".";

    ff:=Conductor(EquationOrder(A));
    T:=MultiplicatorRing(ff);
    assert T in A`KnownOrders;
    // these two attributes are assigned when we created OA above.
    assert assigned T`IsMaximal and assigned T`IsProductOfOrders;
                                                                  // hence T points to the same
                                                                  // memory spot of OA
    printf ".";

    O:=Order(ZBasis(OA));
    assert IsProductOfOrders(O);
    assert IsMaximal(O);
    printf ".";

    O:=MaximalOrder(A);
    G:=[[ Random(O) : i in [1..3] ] : i in [1..100]];
    S:=[ Order(s) : s in G ];
    _:=#Seqset(S);
    printf ".";

    assert forall{z : z in ZBasis(O1) | z in O1 };
    for O in [O1,O2,O3] do
        for i in [1..100] do
            elt:=Random(O);
            assert elt in O;
            abs:=AbsoluteCoordinates(elt,O);
            assert abs eq AbsoluteCoordinates([elt],O)[1];
            assert elt eq DotProduct(abs,ZBasis(O));
        end for;
    end for;
    printf ".";

    seq:=[x^2-5*25,x^2-7*49];
    seq:=[NumberField(f) : f in seq];
    A:=EtaleAlgebra(seq);
    O1:=Order(Basis(A));
    O2:=Order(AbsoluteBasis(A));
    assert O1 eq O2;
    for O in [O1,O2] do
        for i in [1..100] do
            assert Random(O) in O;
        end for;
    end for;
    printf ".";

    E1:=EquationOrder(A);
    E2:=ProductOfEquationOrders(A);
    E3:=Order(A,<EquationOrder(seq[1]) , MaximalOrder(seq[2])>);
    E4:=Order(A,<MaximalOrder(seq[1]) , EquationOrder(seq[2])>);
    OA:=MaximalOrder(A);
    for E in [E1,E2,E3,E4] do
        _:=Index(OA,E);
    end for;
    assert E1 subset E2;
    assert E1 subset E3;
    assert E3*E4 eq OA;
    assert OA meet E1 eq E1;
    assert not E2 meet E1 eq E2;
    printf ".";

    // testing Hash: sets, indexed sets and associative arrays of orders, also when the hash values
    // collide
    _<x>:=PolynomialRing(Integers());
    A:=EtaleAlgebra((x^2+3)*(x^2+27));
    oo:=FindOverOrders(EquationOrder(A)); // distinct orders
    assert #oo eq 85;
    assert #Seqset(oo) eq #oo and #{@ S : S in oo @} eq #oo;
    // copies: the same orders, as new objects whose hash is computed from scratch
    copies1:=[ Order(ZBasis(S) : Check:=0, CheckIsKnownOrder:=false) : S in oo ];
    copies2:=[ Order(ZBasis(S) cat [One(A)] : CheckIsKnownOrder:=false) : S in oo ];
    assert forall{ i : i in [1..#oo] | Hash(copies1[i]) eq Hash(oo[i]) and Hash(copies2[i]) eq Hash(oo[i]) };
    X:={@ S : S in oo @};
    assert forall{ i : i in [1..#oo] | Index(X,copies1[i]) eq i and copies2[i] in X };
    assert #Seqset(oo cat copies1 cat copies2) eq #oo;
    AA:=AssociativeArray();
    for i in [1..#oo] do AA[oo[i]]:=i; end for;
    assert #Keys(AA) eq #oo;
    assert forall{ i : i in [1..#oo] | AA[copies1[i]] eq i and AA[copies2[i]] eq i };
    printf ".";
    // Forced collisions: all the copies get the same hash value, stored in the attribute HashValue.
    // We do not modify the orders in oo, since they are stored in A`KnownOrders.
    for S in copies1 cat copies2 do T:=S; T`HashValue:=0; end for; // T and S are the same object
    assert #{ Hash(S) : S in copies1 cat copies2 } eq 1;
    assert #Seqset(copies1) eq #oo;
    assert #Seqset(copies1 cat copies2) eq #oo;
    X:={@ S : S in copies1 @};
    assert forall{ i : i in [1..#oo] | Index(X,copies2[i]) eq i };
    AA:=AssociativeArray();
    for i in [1..#oo] do AA[copies1[i]]:=i; end for;
    assert #Keys(AA) eq #oo;
    assert forall{ i : i in [1..#oo] | AA[copies2[i]] eq i };
    printf ".";

    // components given by QuadraticField and CyclotomicField: their maximal orders have types
    // RngQuad and RngCyc
    A:=EtaleAlgebra([QuadraticField(-5),CyclotomicField(5)]);
    O:=MaximalOrder(A);
    assert IsMaximal(O);
    test,OasProd:=IsProductOfOrders(O);
    assert test and [ Type(S) : S in OasProd ] eq [ RngQuad, RngCyc ];
    P:=PicardGroup(O);
    assert #P eq 2; // the class numbers of Q(sqrt(-5)) and Q(zeta_5) are 2 and 1
    printf ".";

    // The hashes stored without computing them (from an HNF computed anyway, or copied from an
    // object with the same lattice) must be the ones computed from scratch: we compare with fresh
    // copies, built from the Z-bases, which compute their hash. Intersections of orders and of
    // ideals, ZBasisLLL, restriction of ideals, OneIdeal.
    _<x>:=PolynomialRing(Integers());
    A:=EtaleAlgebra((x^2+5)*(x^4+16));
    E:=EquationOrder(A);
    oo:=OverOrders(E);
    fresh_order:=function(S)
        F:=Order(ZBasis(S) : Check:=0, CheckIsKnownOrder:=false);
        delete F`Hash; // stored by ZBasisLLL at creation: we want it computed from scratch by
                       // myHash
        return F;
    end function;
    fresh_ideal:=func< I | Ideal(Order(I),ZBasis(I)) >; // its hash is computed by ZBasis or myHash
    for i,j in [1..#oo by 3] do
        S:=oo[i] meet oo[j];
        assert myHash(S) eq myHash(fresh_order(S));
    end for;
    pp:=SingularPrimes(E);
    ids:=[ P^2 : P in pp ] cat [ Conductor(E), 2*OneIdeal(E), Ideal(E,[A!7,PrimitiveElement(A)+1]) ];
    for I,J in ids do
        nested:=I subset J or J subset I; // then meet returns one of them, as it is
        K:=I meet J;
        assert (nested or assigned K`Hash) and myHash(K) eq myHash(fresh_ideal(K));
    end for;
    for I in ids do
        J:=fresh_ideal(I);
        // a Z-basis but no hash: ZBasisLLL computes the hash and stores it (a copy of the record of
        // the Z-basis of I, not LLL-reduced)
        Z:=I`ZBasis; Z`lll:=false; J`ZBasis:=Z;
        assert not assigned J`Hash;
        ZBasisLLL(J);
        assert assigned J`Hash and myHash(J) eq myHash(fresh_ideal(J)) and J eq I;
        K:=E!!(oo[#oo]!!I); // the restriction keeps the hash
        assert myHash(K) eq myHash(fresh_ideal(K));
    end for;
    for S in oo do
        assert myHash(OneIdeal(S)) eq myHash(fresh_ideal(OneIdeal(S))) and myHash(OneIdeal(S)) eq myHash(S);
    end for;
    printf ".";

    // Check:=1 suffices when one enlargement is needed; EquationOrder needs an integral primitive
    // element
    _<x>:=PolynomialRing(Integers());
    A:=EtaleAlgebra(x^3-2);
    assert Order([PrimitiveElement(A)] : Check:=1) eq EquationOrder(A);
    _<y>:=PolynomialRing(Rationals());
    failed:=false; try _:=EquationOrder(EtaleAlgebra(y^2-1/2)); catch e failed:=true; end try; assert failed;
    printf ".";
    // with Check:=0, fewer elements than the dimension are refused
    A:=EtaleAlgebra(x^3-2);
    failed:=false; try _:=Order([One(A)] : Check:=0); catch e failed:=true; end try; assert failed;
    printf ".";
    SetAssertions(1);
    printf "all good!";
*/
