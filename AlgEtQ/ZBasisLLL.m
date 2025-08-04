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
// Copyright 2023, S. Marseglia
/////////////////////////////////////////////////////

freeze;

declare verbose ZBasisLLL, 2;

declare attributes AlgEtQIdl : IsZBasisLLLReduced;
declare attributes AlgEtQOrd : IsZBasisLLLReduced;

import "Ord.m" : MatrixAtoQ , MatrixQtoA;

//------------
// LLL - Reduce ZBasis. This should be called whenever we are storing some orders or ideals.
//------------

intrinsic ZBasisLLL(S::AlgEtQOrd)
{A procedure that replaces the ZBasis with an LLL-reduced one. Note: the attribute inclusion matrix, which depends on the Z-Basis is modified as well.}
  if (not assigned S`IsZBasisLLLReduced) or (not S`IsZBasisLLLReduced) then
    S`ZBasis:=MatrixQtoA(Algebra(S),LLL(MatrixAtoQ(ZBasis(S))));
    S`IsZBasisLLLReduced:=true;
    delete S`inclusion_matrix; // this is computed with respect to the old ZBasis!
  end if;
end intrinsic;

intrinsic ZBasisLLL(S::AlgEtQIdl)
{A procedure that replaces the ZBasis with an LLL-reduced one. Note: the attribute inclusion matrix, which depends on the Z-Basis is modified as well.}
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
