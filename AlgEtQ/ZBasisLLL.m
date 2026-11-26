/////////////////////////////////////////////////////
// Copyright 2025.
// Stefano Marseglia, stefano.marseglia89@gmail.com
// https://stmar89.github.io/index.html
//
// Distributed under the terms of the CC-BY 4.0 licence.
// https://creativecommons.org/licenses/by/4.0/
/////////////////////////////////////////////////////



freeze;

declare verbose ZBasisLLL, 2;

import "Ord.m" : MatrixAtoQ , MatrixQtoA, crQZ, crZQ, hash_from_hnf, hnf_from_hash, lattice_hnf,
                 zb_from_seq, zb_from_dM, zb_seq, zb_dM, zb_inv, replace_zbasis_if_smaller;

///# LLL-reduction of $\mathbb{Z}$-bases
/// Whenever we call the intrinsic `ZBasis` for an order or a fractional ideal, we are computing a $\mathbb{Z}$-basis. This is produced directly from the known generators. This elements might have large coefficients making certain operation considerable slower. One can then reduce the given basis using the LLL algorithm to get better coefficients. While this procedure might cost some time, it is often convenient to perform this reduction if the ideal or the order are supposed to be used further.
/// We provide the intrinsic `ZBasisLLL` which is a procedure changing the stored $\mathbb{Z}$-basis to an LLL-reduced one, if this has not been done before.
/// The running time of the LLL algorithm grows with the size of the coefficients of the input basis. If the Hermite normal form of the basis has smaller coefficients, which happens for example when the stored basis has very large coefficients but the lattice has small index, then the LLL algorithm is applied to the Hermite normal form instead.
/// The output of `ZBasisLLL` depends on the stored $\mathbb{Z}$-basis. The intrinsic `ZBasisLLLFromHNF` returns instead the LLL reduction of the Hermite normal form, which depends only on the lattice.
/// The stored $\mathbb{Z}$-basis is also replaced when a smaller $\mathbb{Z}$-basis of the same lattice becomes known, for example the Hermite normal form computed for the hash: see `ReplaceZBasis`.

//------------
// LLL - Reduce ZBasis. This should be called whenever we are storing some orders or ideals.
//------------

max_abs:=function(X)
// the maximum of the absolute values of the entries of the integer matrix X
    E:=Eltseq(X);
    return Max(Max(E),-Min(E));
end function;

lll_zbasis:=function(X)
// Given a lattice X (an order, a fractional ideal or a module) with a working Z-basis (see
// zbasis_fmt in Ord.m), returns the record of an LLL reduction of it. The hash of X is stored if it
// was not known: the HNF is computed anyway.
// Let d and M be as in zbasis_fmt (M = d times the matrix of the coordinates of the basis).
// The entries of the HNF H of M are bounded by its pivots, while the entries of M are not:
// if H has smaller entries, we LLL-reduce it instead, which is much faster.
// Otherwise, for example if the basis is already (almost) reduced, we LLL-reduce M itself,
// since starting from the HNF would undo the reduction already present in the basis.
// Working with the integer matrix M is also faster than calling LLL on a rational matrix.
// H can be read from the hash, if known: otherwise we compute it, and the hash comes for free.
    Z:=X`ZBasis;
    d,M:=zb_dM(X);
    if assigned X`Hash then
        H:=hnf_from_hash(X`Hash,d);
        assert2 H eq HermiteForm(M);
    elif Z`hnf then
        H:=M;
        X`Hash:=hash_from_hnf(H,d);
    else
        H:=HermiteForm(M); // M is square of full rank, so H has no zero rows
        assert2 Rank(H) eq Nrows(H);
        X`Hash:=hash_from_hnf(H,d);
    end if;
    use_hnf:=max_abs(H) lt max_abs(M);
    vprintf ZBasisLLL, 2 : "lll_zbasis: max entry of the basis has %o bits, of the HNF %o bits. Using the %o.\n",
                            Ilog2(max_abs(M)+1), Ilog2(max_abs(H)+1), use_hnf select "HNF" else "basis";
    L:=LLL(use_hnf select H else M);
    assert2 HermiteForm(L) eq H; // same lattice
    return zb_from_dM(Z`A,d,L : lll:=true);
end function;

canonical_zbasis:=function(X)
// Given a lattice X (an order, a fractional ideal or a module), returns the record (see zbasis_fmt
// in Ord.m) of the LLL reduction of its HNF, read from the hash: it depends only on the lattice.
    h:=myHash(X); // also makes sure that X has a working Z-basis
    d:=h[1];
    return zb_from_dM(X`ZBasis`A,d,LLL(hnf_from_hash(h,d)) : lll:=true);
end function;

/// A procedure that replaces the stored $\mathbb{Z}$-basis with an LLL-reduced one, if this has not been done before.
intrinsic ZBasisLLL(S::AlgEtQOrd)
{A procedure that replaces the stored ZBasis with an LLL-reduced one, if this has not been done before.}
  if not S`ZBasis`lll then
    _:=replace_zbasis_if_smaller(S,lll_zbasis(S) : Force:=true);
  end if;
end intrinsic;

///ditto
intrinsic ZBasisLLL(S::AlgEtQIdl)
{ " } // "
  if not assigned S`ZBasis then
    _:=ZBasis(S);
  end if;
  if not S`ZBasis`lll then
    // If the ideal has at least dim generators, for example the old ZBasis when it was created from
    // more than dim generators, they are replaced with the LLL-reduced ZBasis (see
    // replace_zbasis_if_smaller): it generates the ideal, it is not longer, and its elements are
    // smaller, which makes the operations using Generators faster, e.g. products of ideals. Fewer
    // generators are kept, since products of ideals with few generators avoid computing a ZBasis.
    _:=replace_zbasis_if_smaller(S,lll_zbasis(S) : Force:=true);
  end if;
end intrinsic;

/// Returns the LLL reduction of the Hermite normal form of the order (or of the fractional ideal): a $\mathbb{Z}$-basis which depends only on the lattice, and not on the stored $\mathbb{Z}$-basis. It is computed once and stored. If it is smaller than the stored $\mathbb{Z}$-basis, it replaces it, see `ReplaceZBasis`.
intrinsic ZBasisLLLFromHNF(S::AlgEtQOrd) -> SeqEnum[AlgEtQElt]
{Returns the LLL reduction of the Hermite normal form of the order (or of the fractional ideal): a Z-basis which depends only on the lattice, and not on the stored Z-basis. It is computed once and stored. If it is smaller than the stored Z-basis, it replaces it, see ReplaceZBasis.}
    if not assigned S`CanonicalZBasis then
        S`CanonicalZBasis:=canonical_zbasis(S);
    end if;
    zb:=zb_seq(S : att:="CanonicalZBasis"); // before the record is copied below
    _:=replace_zbasis_if_smaller(S,S`CanonicalZBasis);
    return zb;
end intrinsic;

///ditto
intrinsic ZBasisLLLFromHNF(S::AlgEtQIdl) -> SeqEnum[AlgEtQElt]
{ " } // "
    if not assigned S`CanonicalZBasis then
        S`CanonicalZBasis:=canonical_zbasis(S);
    end if;
    zb:=zb_seq(S : att:="CanonicalZBasis"); // before the record is copied below
    _:=replace_zbasis_if_smaller(S,S`CanonicalZBasis);
    return zb;
end intrinsic;

/// Replaces the stored $\mathbb{Z}$-basis of $S$ with the given $\mathbb{Z}$-basis `zb` of the same lattice, if it is smaller: if the integer matrix of its coordinates (multiplied by their common denominator) has at most $3/4$ of the bits of the one of the stored $\mathbb{Z}$-basis. With `Force:=true` it is always replaced. The stored $\mathbb{Z}$-basis of orders, ideals and modules is replaced only in this way (also by the package itself, e.g. when the Hermite normal form is computed for the hash, and by `ZBasisLLL`): with `SetVerbose("ZBasis",1)` every replacement which is not forced is printed.
intrinsic ReplaceZBasis(S::AlgEtQOrd, zb::SeqEnum[AlgEtQElt] : Force:=false)
{Replaces the stored Z-basis of S with the given Z-basis zb of the same lattice, if it is smaller: if the integer matrix of its coordinates (multiplied by their common denominator) has at most 3/4 of the bits of the one of the stored Z-basis. With Force:=true it is always replaced. The stored Z-basis of orders, ideals and modules is replaced only in this way (also by the package itself, e.g. when the Hermite normal form is computed for the hash, and by ZBasisLLL): with SetVerbose("ZBasis",1) every replacement which is not forced is printed.}
    A:=Algebra(S);
    require #zb eq AbsoluteDimension(A) and Universe(zb) cmpeq A : "zb must be a Z-basis in the algebra of S";
    H,d:=lattice_hnf(zb);
    require Nrows(H) eq #zb and hash_from_hnf(H,d) eq myHash(S) : "zb must be a Z-basis of S";
    _:=replace_zbasis_if_smaller(S,zb_from_seq(A,zb) : Force:=Force);
end intrinsic;

///ditto
intrinsic ReplaceZBasis(S::AlgEtQIdl, zb::SeqEnum[AlgEtQElt] : Force:=false)
{ " } // "
    A:=Algebra(S);
    require #zb eq AbsoluteDimension(A) and Universe(zb) cmpeq A : "zb must be a Z-basis in the algebra of S";
    H,d:=lattice_hnf(zb);
    require Nrows(H) eq #zb and hash_from_hnf(H,d) eq myHash(S) : "zb must be a Z-basis of S";
    _:=replace_zbasis_if_smaller(S,zb_from_seq(A,zb) : Force:=Force);
end intrinsic;

//------------
// Short elements modulo a fractional ideal
//------------

/// Given an element $x$ and a fractional ideal $J$ in the same algebra, returns an element of $x+J$ with small coordinates. It is obtained by subtracting from $x$ the element of $J$ whose coordinates with respect to the stored $\mathbb{Z}$-basis of $J$ are the coordinates of $x$ rounded to the nearest integers. The stored $\mathbb{Z}$-basis of $J$ is first replaced by an LLL-reduced one, see `ZBasisLLL`. With the vararg `Canonical:=true` (default `false`) the $\mathbb{Z}$-basis `ZBasisLLLFromHNF(J)` is used instead, which depends only on $J$: then the output depends only on $x+J$, and not on the stored $\mathbb{Z}$-basis of $J$. In both cases the output depends only on the class of $x$ modulo $J$ (given the $\mathbb{Z}$-basis), and it is $0$ if and only if $x$ is in $J$.
intrinsic ShortElementModulo(x::AlgEtQElt, J::AlgEtQIdl : Canonical:=false) -> AlgEtQElt
{Given an element x and a fractional ideal J in the same algebra, returns an element of x+J with small coordinates. It is obtained by subtracting from x the element of J whose coordinates with respect to the stored Z-basis of J are the coordinates of x rounded to the nearest integers. The stored Z-basis of J is first replaced by an LLL-reduced one, see ZBasisLLL. With the vararg Canonical:=true (default false) the Z-basis ZBasisLLLFromHNF(J) is used instead, which depends only on J: then the output depends only on x+J, and not on the stored Z-basis of J. In both cases the output depends only on the class of x modulo J (given the Z-basis), and it is 0 if and only if x is in J.}
    require Algebra(x) cmpeq Algebra(J) : "the element and the ideal must be in the same algebra";
    if Canonical then
        zb:=ZBasisLLLFromHNF(J);
        Minv:=zb_inv(J : att:="CanonicalZBasis");
    else
        ZBasisLLL(J);
        zb:=ZBasis(J);
        Minv:=zb_inv(J);
    end if;
    c:=Eltseq(MatrixAtoQ([x])*Minv); // the coordinates of x with respect to zb
    // We round with Floor(e+1/2), since it commutes with adding integers: hence the output depends
    // only on x+J.
    return x - DotProduct([ Floor(e+1/2) : e in c ], zb);
end intrinsic;

/* TESTS

  printf "### Testing ZBasisLLL:";
    //AttachSpec("../spec");
    SetAssertions(2);
    _<x>:=PolynomialRing(Integers());
  f:=x^4-100*x^3-100*x^2-100*x-100;
  K:=EtaleAlgebra(f);
  E:=EquationOrder(K);
  pp:=SingularPrimes(E);
  I:=&*(pp);
  J:=&*(pp);
  ZBasisLLL(I);
  assert ZBasis(J) ne ZBasis(I);
  assert J eq I;

  // helpers
  zbl_AtoQ:=func< gens | Matrix(Rationals(),[ AbsoluteCoordinates(g) : g in gens ]) >;
  zbl_QtoA:=func< A,M | [ &+[ M[r,i]*AbsoluteBasis(A)[i] : i in [1..Ncols(M)] ] : r in [1..Nrows(M)] ] >;
  zbl_max:=func< X | Max([ Abs(c) : c in Eltseq(X) ]) >;
  zbl_is_lll_reduced:=function(M)
  // exact Gram-Schmidt over Q: size reduction with eta=0.51 and Lovasz condition with delta=3/4
      n:=Nrows(M);
      Bs:=[]; nn:=[]; mu:=ZeroMatrix(Rationals(),n,n);
      for i in [1..n] do
          v:=M[i];
          for j in [1..i-1] do
              mu[i,j]:=InnerProduct(M[i],Bs[j])/nn[j];
              v:=v-mu[i,j]*Bs[j];
          end for;
          Append(~Bs,v); Append(~nn,InnerProduct(v,v));
      end for;
      return forall{ 0 : j in [1..i-1], i in [1..n] | Abs(mu[i,j]) le 51/100 } and
             forall{ i : i in [2..n] | nn[i] ge (3/4-mu[i,i-1]^2)*nn[i-1] };
  end function;
  zbl_skew:=function(M,b)
  // left multiplication by a random unimodular matrix with large entries
      n:=Nrows(M); U:=IdentityMatrix(Integers(),n);
      for t in [1..2*n] do
          i:=Random(1,n); repeat j:=Random(1,n); until j ne i;
          AddRow(~U,Random(-2^b,2^b),j,i);
      end for;
      return ChangeRing(U,Rationals())*M;
  end function;
  zbl_check:=procedure(S,zb,branch)
  // S is an order or an ideal with stored ZBasis zb, not flagged as LLL-reduced.
  // branch = "HNF", "basis" or "any": the matrix that ZBasisLLL should LLL-reduce.
      A:=Algebra(S);
      M:=zbl_AtoQ(zb);
      d:=Denominator(M);
      MZ:=ChangeRing(d*M,Integers());
      H:=HermiteForm(MZ);
      use_hnf:=zbl_max(H) lt zbl_max(MZ);
      if branch eq "HNF" then assert use_hnf; elif branch eq "basis" then assert not use_hnf; end if;
      ZBasisLLL(S);
      MS:=zbl_AtoQ(ZBasis(S));
      // the expected branch was used
      assert MS eq (1/d)*ChangeRing(LLL(use_hnf select H else MZ),Rationals());
      assert HermiteForm(ChangeRing(d*MS,Integers())) eq H; // same lattice
      assert zbl_is_lll_reduced(MS);
      assert S`ZBasis`lll;
      zb1:=ZBasis(S);
      ZBasisLLL(S);
      assert ZBasis(S) eq zb1; // nothing changes the second time
      // membership uses the inclusion matrix, which must be the one of the new ZBasis
      assert forall{ z : z in zb | z in S };
      assert forall{ i : i in [1..3] | DotProduct([ Random(-5,5) : z in zb1 ],zb1) in S };
      if Type(S) eq AlgEtQOrd then
          assert One(A) in S;
      end if;
  end procedure;
  zbl_fresh:=function(R,zb) J:=Ideal(R,zb); ReplaceZBasis(J,zb : Force:=true); return J; end function;

  SetSeed(1);
  pols:=[ x^8+16, (x-1)*(x-20)*(x^2-5)*(x^3-2), x^10-x^9+4*x^8-6*x^7+8*x^6-16*x^5+16*x^4-24*x^3+32*x^2-16*x+32 ];
  for f in pols do
      A:=EtaleAlgebra(f);
      R:=EquationOrder(A);
      N:=AbsoluteDimension(A);
      repeat y:=A!Random(R,3); until not IsZeroDivisor(y);
      I:=Ideal(R,[ (A!Random(R,3))/y : i in [1..3] ]);
      I:=I*Ideal(R,[ A!Random(R,3) : i in [1..3] ]);
      zbI:=ZBasis(I);
      // 1. HNF basis
      zbl_check(zbl_fresh(R,zbI),zbI,"any");
      printf ".";
      // 2. skewed basis: huge entries, small index
      zb:=zbl_QtoA(A,zbl_skew(zbl_AtoQ(zbI),40));
      zbl_check(zbl_fresh(R,zb),zb,"HNF");
      printf ".";
      // 3. already LLL-reduced basis, not flagged
      zb:=zbl_QtoA(A,LLL(zbl_AtoQ(zbI)));
      zbl_check(zbl_fresh(R,zb),zb,"basis");
      printf ".";
      // 4. huge index, and its colon
      Ik:=I^(N le 8 select 20 else 8);
      zbl_check(zbl_fresh(R,ZBasis(Ik)),ZBasis(Ik),"any");
      C:=ColonIdeal(R,Ik);
      zbl_check(zbl_fresh(R,ZBasis(C)),ZBasis(C),"any");
      printf ".";
      // 5. a*I as in SmallRepresentative
      a:=ShortElement(ColonIdeal(R,Ik));
      zb:=[ a*z : z in ZBasis(Ik) ];
      zbl_check(zbl_fresh(R,zb),zb,"HNF");
      printf ".";
      // orders: the maximal order with a skewed basis
      O:=MaximalOrder(A);
      zb:=zbl_QtoA(A,zbl_skew(zbl_AtoQ(ZBasis(O)),40));
      ReplaceZBasis(O,zb : Force:=true); // a record which is not flagged as LLL-reduced
      zbl_check(O,zb,"HNF");
      assert R subset O;
      printf ".";
  end for;

  // ZBasisLLLFromHNF: the LLL reduction of the HNF, which depends only on the lattice and not on
  // the stored Z-basis; it replaces the stored Z-basis only if it is smaller (see ReplaceZBasis).
  A:=EtaleAlgebra(x^8+16);
  R:=EquationOrder(A);
  I:=Ideal(R,[ A!Random(R,3) : i in [1..3] ])^20;
  zbI:=ZBasis(I); // the HNF basis
  I1:=zbl_fresh(R,zbI);
  c1:=ZBasisLLLFromHNF(I1);
  ZBasisLLL(I1);
  assert ZBasis(I1) eq c1; // the LLL reduction of the HNF basis is the canonical one
  zb:=zbl_QtoA(A,LLL(zbl_skew(zbl_AtoQ(zbI),5))); // reduced, from another basis
  I2:=zbl_fresh(R,zb);
  ZBasisLLL(I2);
  zb2:=ZBasis(I2);
  ZBasisLLL(I2);
  assert ZBasis(I2) eq zb2; // already reduced: nothing to do
  assert ZBasisLLLFromHNF(I2) eq c1; // the same lattice
  I3:=zbl_fresh(R,zbl_QtoA(A,zbl_skew(zbl_AtoQ(zbI),40))); // huge entries
  assert ZBasisLLLFromHNF(I3) eq c1;
  assert ZBasis(I3) eq c1; // much smaller than the stored Z-basis: it replaced it
  assert forall{ z : z in zbI | z in I2 and z in I3 }; // the inclusion matrices are the new ones
  // orders: the canonical Z-basis of the maximal order, and orders with other Z-bases
  O:=MaximalOrder(A);
  cO:=ZBasisLLLFromHNF(O);
  // Order reduces the HNF
  O2:=Order(zbl_QtoA(A,zbl_skew(zbl_AtoQ(cO),5)) : CheckIsKnownOrder:=false);
  assert ZBasis(O2) eq cO and ZBasisLLLFromHNF(O2) eq cO and One(A) in O2;
  ReplaceZBasis(O2,zbl_QtoA(A,zbl_skew(zbl_AtoQ(cO),40)) : Force:=true);
  // smaller: it replaced it
  assert ZBasisLLLFromHNF(O2) eq cO and ZBasis(O2) eq cO and One(A) in O2;
  // without Force, a Z-basis which is not smaller is not taken
  zbO2:=ZBasis(O2);
  ReplaceZBasis(O2,zbl_QtoA(A,zbl_skew(zbl_AtoQ(cO),40)));
  assert ZBasis(O2) eq zbO2;
  // the HNF computed for the hash replaces a Z-basis with many more bits
  J:=zbl_fresh(R,zbl_QtoA(A,zbl_skew(zbl_AtoQ(zbI),40)));
  delete J`Hash; if assigned J`HashValue then delete J`HashValue; end if; // as if not known
  bJ:=ZBasis(J);
  _:=myHash(J);
  assert ZBasis(J) eq zbI and J eq Ideal(R,bJ); // zbI is the HNF basis
  printf ".";

  // ShortElementModulo(y,J) is in y+J, it depends only on y+J, it is 0 on J, and its coordinates
  // are at most half the sum of the absolute values of the coordinates of the elements of the
  // LLL-reduced Z-basis of J.
  A:=EtaleAlgebra(x^8+16);
  E:=EquationOrder(A);
  pp:=SingularPrimes(E);
  ids:=[ OneIdeal(E), Conductor(E), 3*OneIdeal(E), (1/2)*Conductor(E) ] cat pp cat [ P^3 : P in pp ];
  for J in ids do
      for i in [1..5] do
          y:=DotProduct([ Random(-10^10,10^10) : j in [1..8] ],AbsoluteBasis(A))/Random(1,100);
          z:=ShortElementModulo(y,J);
          assert y-z in J;
          assert ShortElementModulo(y+Random(J,10^5),J) eq z;
          bound:=[ (1/2)*&+[ Abs(AbsoluteCoordinates(b)[k]) : b in ZBasis(J) ] : k in [1..8] ];
          assert forall{ k : k in [1..8] | Abs(AbsoluteCoordinates(z)[k]) le bound[k] };
      end for;
      assert ShortElementModulo(Random(J,100),J) eq Zero(A);
      printf ".";
  end for;

  // ZBasisLLL replaces the generators of an ideal with the LLL-reduced ZBasis when they are at
  // least dim, and keeps them when they are fewer.
  A:=EtaleAlgebra(x^8+16);
  E:=EquationOrder(A);
  for k in [2,8,12] do
      I:=Ideal(E,[ Random(E,10^6) : i in [1..k] ]);
      gens:=Generators(I);
      J:=Ideal(E,gens);
      ZBasisLLL(I);
      assert I eq J;
      if k ge 8 then
          assert Generators(I) eq ZBasis(I);
      else
          assert Generators(I) eq gens;
      end if;
      assert Ideal(E,Generators(I)) eq J;
      assert I*I eq J*J and I+OneIdeal(E) eq J+OneIdeal(E);
  end for;
  printf ".";

  SetAssertions(1);
  printf " all good!";

*/
