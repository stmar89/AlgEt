/////////////////////////////////////////////////////
// Stefano Marseglia, stefano.marseglia89@gmail.com
// https://stmar89.github.io/index.html
// 
// Distributed under the terms of the GNU Lesser General Public License (L-GPL)
//      http://www.gnu.org/licenses/
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 3.0 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301  USA
// 
// Copyright 2025, Stefano Marseglia
/////////////////////////////////////////////////////


freeze;

declare verbose ZBasisLLL, 2;

declare attributes AlgEtQIdl : IsZBasisLLLReduced;
declare attributes AlgEtQOrd : IsZBasisLLLReduced;

import "Ord.m" : MatrixAtoQ , MatrixQtoA;

///# LLL-reduction of $\mathbb{Z}$-bases
/// Whenever we call the intrinsic `ZBasis` for an order or a fractional ideal, we are computing a $\mathbb{Z}$-basis. This is produced directly from the known generators. This elements might have large coefficients making certain operation considerable slower. One can then reduce the given basis using the LLL algorithm to get better coefficients. While this procedure might cost some time, it is often convenient to perform this reduction if the ideal or the order are supposed to be used further.
/// We provide the intrinsic `ZBasisLLL` which is a procedure changing the stored $\mathbb{Z}$-basis to an LLL-reduced one, if this has not been done before. The attribute `inclusion_matrix` is also deleted, since it depends on the stored $\mathbb{Z}$-basis.

//------------
// LLL - Reduce ZBasis. This should be called whenever we are storing some orders or ideals.
//------------

/// A procedure that replaces the stored $\mathbb{Z}$-basis with an LLL-reduced one. 
intrinsic ZBasisLLL(S::AlgEtQOrd)
{A procedure that replaces the stored ZBasis with an LLL-reduced one. Note: the attribute inclusion matrix, which depends on the ZBasis is deleted.}
  if (not assigned S`IsZBasisLLLReduced) or (not S`IsZBasisLLLReduced) then
    S`ZBasis:=MatrixQtoA(Algebra(S),LLL(MatrixAtoQ(ZBasis(S))));
    S`IsZBasisLLLReduced:=true;
    delete S`inclusion_matrix; // this is computed with respect to the old ZBasis!
  end if;
end intrinsic;

///ditto
intrinsic ZBasisLLL(S::AlgEtQIdl)
{A procedure that replaces the stored ZBasis with an LLL-reduced one. Note: the attribute inclusion matrix, which depends on the ZBasis is deleted.}
  if (not assigned S`IsZBasisLLLReduced) or (not S`IsZBasisLLLReduced) then
    S`ZBasis:=MatrixQtoA(Algebra(S),LLL(MatrixAtoQ(ZBasis(S))));
    S`IsZBasisLLLReduced:=true;
    delete S`inclusion_matrix; // this is computed with respect to the old ZBasis!
  end if;
end intrinsic;

/* TESTS

  printf "### Testing ZBasisLLL:";
	//AttachSpec("~/packages_github/AlgEt/spec");
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
  printf " all good!"; 
*/