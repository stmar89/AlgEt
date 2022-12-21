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
  if not assigned S`IsZBasisLLLReduced or not S`IsZBasisLLLReduced then
    S`ZBasis:=MatrixQtoA(Algebra(S),LLL(MatrixAtoQ(ZBasis(S))));
    S`IsZBasisLLLReduced:=true;
  end if;
end intrinsic;

intrinsic ZBasisLLL(S::AlgEtIdl)
{A procedure that replaces the ZBasis with an LLL-reduced one.}
  if not assigned S`IsZBasisLLLReduced or not S`IsZBasisLLLReduced then
    S`ZBasis:=MatrixQtoA(Algebra(S),LLL(MatrixAtoQ(ZBasis(S))));
    S`IsZBasisLLLReduced:=true;
  end if;
end intrinsic;



/* TEST

*/
