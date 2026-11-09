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

declare attributes AlgEtQIdl : IsZBasisLLLReduced;
declare attributes AlgEtQOrd : IsZBasisLLLReduced;

import "Ord.m" : MatrixAtoQ , MatrixQtoA, crQZ, crZQ;

///# LLL-reduction of $\mathbb{Z}$-bases
/// Whenever we call the intrinsic `ZBasis` for an order or a fractional ideal, we are computing a $\mathbb{Z}$-basis. This is produced directly from the known generators. This elements might have large coefficients making certain operation considerable slower. One can then reduce the given basis using the LLL algorithm to get better coefficients. While this procedure might cost some time, it is often convenient to perform this reduction if the ideal or the order are supposed to be used further.
/// We provide the intrinsic `ZBasisLLL` which is a procedure changing the stored $\mathbb{Z}$-basis to an LLL-reduced one, if this has not been done before. The attribute `inclusion_matrix` is also deleted, since it depends on the stored $\mathbb{Z}$-basis.
/// The running time of the LLL algorithm grows with the size of the coefficients of the input basis. If the Hermite normal form of the basis has smaller coefficients, which happens for example when the stored basis has very large coefficients but the lattice has small index, then the LLL algorithm is applied to the Hermite normal form instead.

//------------
// LLL - Reduce ZBasis. This should be called whenever we are storing some orders or ideals.
//------------

max_abs:=function(X)
// the maximum of the absolute values of the entries of the integer matrix X
    E:=Eltseq(X);
    return Max(Max(E),-Min(E));
end function;

lll_zbasis:=function(A,zb)
// Given a Z-basis zb of a lattice in A, returns an LLL-reduced Z-basis of the same lattice.
// Let M be the matrix of the coordinates of zb and d its denominator.
// The entries of the HNF of d*M are bounded by its pivots, while the entries of d*M are not:
// if the HNF has smaller entries, we LLL-reduce it instead, which is much faster.
// Otherwise, for example if zb is already (almost) reduced, we LLL-reduce d*M itself,
// since starting from the HNF would undo the reduction already present in zb.
// Working with the integer matrix d*M is also faster than calling LLL on M.
    M:=MatrixAtoQ(zb);
    d:=Denominator(M);
    MZ:=crQZ(d*M);
    H:=HermiteForm(MZ); // MZ is square of full rank, so H has no zero rows
    assert2 Rank(H) eq Nrows(H);
    use_hnf:=max_abs(H) lt max_abs(MZ);
    vprintf ZBasisLLL, 2 : "lll_zbasis: max entry of the basis has %o bits, of the HNF %o bits. Using the %o.\n",
                            Ilog2(max_abs(MZ)+1), Ilog2(max_abs(H)+1), use_hnf select "HNF" else "basis";
    if use_hnf then
        MZ:=H;
    end if;
    L:=LLL(MZ);
    assert2 HermiteForm(L) eq H; // same lattice
    return MatrixQtoA(A,(1/d)*crZQ(L));
end function;

/// A procedure that replaces the stored $\mathbb{Z}$-basis with an LLL-reduced one.
intrinsic ZBasisLLL(S::AlgEtQOrd)
{A procedure that replaces the stored ZBasis with an LLL-reduced one. Note: the attribute inclusion matrix, which depends on the ZBasis is deleted.}
  if (not assigned S`IsZBasisLLLReduced) or (not S`IsZBasisLLLReduced) then
    S`ZBasis:=lll_zbasis(Algebra(S),ZBasis(S));
    S`IsZBasisLLLReduced:=true;
    delete S`inclusion_matrix; // this is computed with respect to the old ZBasis!
  end if;
end intrinsic;

///ditto
intrinsic ZBasisLLL(S::AlgEtQIdl)
{ " } // "
  if (not assigned S`IsZBasisLLLReduced) or (not S`IsZBasisLLLReduced) then
    S`ZBasis:=lll_zbasis(Algebra(S),ZBasis(S));
    S`IsZBasisLLLReduced:=true;
    delete S`inclusion_matrix; // this is computed with respect to the old ZBasis!
  end if;
end intrinsic;

//------------
// Short elements modulo a fractional ideal
//------------

/// Given an element $x$ and a fractional ideal $J$ in the same algebra, returns an element of $x+J$ with small coordinates. It is obtained by subtracting from $x$ the element of $J$ whose coordinates with respect to the stored $\mathbb{Z}$-basis of $J$ are the coordinates of $x$ rounded to the nearest integers. The stored $\mathbb{Z}$-basis of $J$ is first replaced by an LLL-reduced one, see `ZBasisLLL`. The output depends only on the class of $x$ modulo $J$, and it is $0$ if and only if $x$ is in $J$.
intrinsic ShortElementModulo(x::AlgEtQElt, J::AlgEtQIdl) -> AlgEtQElt
{Given an element x and a fractional ideal J in the same algebra, returns an element of x+J with small coordinates. It is obtained by subtracting from x the element of J whose coordinates with respect to the stored Z-basis of J are the coordinates of x rounded to the nearest integers. The stored Z-basis of J is first replaced by an LLL-reduced one, see ZBasisLLL. The output depends only on the class of x modulo J, and it is 0 if and only if x is in J.}
    require Algebra(x) cmpeq Algebra(J) : "the element and the ideal must be in the same algebra";
    ZBasisLLL(J);
    c := AbsoluteCoordinates([x],J)[1];
    // We round with Floor(e+1/2), since it commutes with adding integers: hence the output depends only on x+J.
    return x - DotProduct([ Floor(e+1/2) : e in c ], ZBasis(J));
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
      assert MS eq (1/d)*ChangeRing(LLL(use_hnf select H else MZ),Rationals()); // the expected branch was used
      assert HermiteForm(ChangeRing(d*MS,Integers())) eq H; // same lattice
      assert zbl_is_lll_reduced(MS);
      assert S`IsZBasisLLLReduced;
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
  zbl_fresh:=function(R,zb) J:=Ideal(R,zb); J`ZBasis:=zb; return J; end function;

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
      O`ZBasis:=zb;
      O`IsZBasisLLLReduced:=false;
      delete O`inclusion_matrix;
      zbl_check(O,zb,"HNF");
      assert R subset O;
      printf ".";
  end for;

  // ShortElementModulo(y,J) is in y+J, it depends only on y+J, it is 0 on J, and its coordinates are at most half
  // the sum of the absolute values of the coordinates of the elements of the LLL-reduced Z-basis of J.
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

  SetAssertions(1);
  printf " all good!";

*/
