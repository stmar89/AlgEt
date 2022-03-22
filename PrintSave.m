/* vim: set syntax=magma :*/

freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

//declare verbose ?????, 1;

/*TODO:

*/

//------------
// Print to string
//------------

intrinsic PrintSeqAlgEtElt(seq::SeqEnum[AlgEtElt]) -> SeqEnum
{ Given a sequence of elements of an AlgEt, returns a sequence of tuples of sequence of integers that can be coerced into the original algebra to obtain the input sequece.
  Such an output can easily be printed to file.}
  return [ < Eltseq(c) : c in Components(elt) > : elt in seq ];
end intrinsic;




/* TEST
    
    AttachSpec("~/packages_github/AlgEt/spec");
    _<x>:=PolynomialRing(Integers());
    f:=x^6 - 3*x^5 - 3*x^4 + 65*x^3 - 48*x^2 - 768*x + 4096;
    A:=EtaleAlgebra(f);
    E:=EquationOrder(A);
    seq:=PrintSeqAlgEtElt(ZBasis(E));
    seq:=[ A! s : s in seq ];
    assert Order(seq) eq E;




*/
