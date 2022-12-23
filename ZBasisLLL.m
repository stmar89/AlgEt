/* vim: set syntax=magma :*/

freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose ZBasisLLL, 2;

declare attributes AlgEtIdl : IsZBasisLLLReduced;
declare attributes AlgEtOrd : IsZBasisLLLReduced;

import "Ord.m" : MatrixAtoQ , MatrixQtoA;

//------------
// LLL - Reduce ZBasis. This should be called whenever we are storing some orders or ideals.
//------------

intrinsic ZBasisLLL(S::AlgEtOrd)
{A procedure that replaces the ZBasis with an LLL-reduced one.}
  if (not assigned S`IsZBasisLLLReduced) or (not S`IsZBasisLLLReduced) then
    S`ZBasis:=MatrixQtoA(Algebra(S),LLL(MatrixAtoQ(ZBasis(S))));
    S`IsZBasisLLLReduced:=true;
    delete S`inclusion_matrix; // this is computed with respect to the old ZBasis!
  end if;
end intrinsic;

intrinsic ZBasisLLL(S::AlgEtIdl)
{A procedure that replaces the ZBasis with an LLL-reduced one.}
  if (not assigned S`IsZBasisLLLReduced) or (not S`IsZBasisLLLReduced) then
    S`ZBasis:=MatrixQtoA(Algebra(S),LLL(MatrixAtoQ(ZBasis(S))));
    S`IsZBasisLLLReduced:=true;
    delete S`inclusion_matrix; // this is computed with respect to the old ZBasis!
  end if;
end intrinsic;

/* TEST

  printf "### Testing ZBasisLLL:";
	AttachSpec("~/packages_github/AlgEt/spec");
	SetAssertions(2);
	_<x>:=PolynomialRing(Integers());
  f:=x^4-100*x^3-100*x^2-100*x-100;
  K:=EtaleAlgebra(f);
  E:=EquationOrder(K);
  pp:=PrimesAbove(Conductor(E));
  I:=&*(pp);
  J:=&*(pp);
  ZBasisLLL(I);
  assert ZBasis(J) ne ZBasis(I);
  assert J eq I;
  printf " all good!\n"; 

*/
