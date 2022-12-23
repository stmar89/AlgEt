/* vim: set syntax=magma :*/

//freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose AlgEtQ, 1;

/*TODO:

*/

declare type AlgEtQ[AlgEtQElt];

declare attributes AlgEtQ : DefiningPolynomial, 
                           // ass_algebra, 
                           Dimension,
                           AbsoluteDimension,
                           BaseField, //a tup : <F,m> where F is the Base field and m is the diagonal embedding into A
                           HasBaseField, //a boolean
                           PrimeField,
                           Components; //a tup of 3 sequences: the first are the NF, 
                                         //the second are embeddings and the third are projections

//------------
// Creation for AlgEtQ
//------------

intrinsic EtaleAlgebra(seq::SeqEnum[FldNum[FldRat]]) -> AlgEtQ
{Given a sequence of number fields returns the étale algebra corresponding to the direct product.}
    A:=New(AlgEtQ);
    embs:=[ map< seq[i]->A | x:-> A! (<seq[j]!0 : j in [1..i-1]> cat <x> cat <seq[j]!0 : j in [i+1..#seq]>)  > : i in [1..#seq] ];
    projs:=[ map< A->seq[i] | y:-> Components(y)[i] > : i in [1..#seq] ];
    A`Components:=<seq,embs,projs>;
    return A;
end intrinsic;

intrinsic EtaleAlgebra(f::RngUPolElt[RngInt]) -> AlgEtQ
{Given a squarefree polynomial over the integers returns the product of the number fields defined by the irreducible factors.}
    require IsSquarefree(f) : "The polynomial must be squarefree.";
    A:=EtaleAlgebra([NumberField(g[1]) : g in Factorization(f)]);
    A`DefiningPolynomial:=f;
    return A;
end intrinsic;

intrinsic EtaleAlgebra(f::RngUPolElt[FldRat]) -> AlgEtQ
{Given a squarefree polynomial over the rationals returns the product of the number fields defined by the irreducible factors.}
    require IsSquarefree(f) : "The polynomial must be squarefree.";
    A:=EtaleAlgebra([NumberField(g[1]) : g in Factorization(f)]);
    A`DefiningPolynomial:=f;
    return A;
end intrinsic;

/* TESTS

    printf "### Testing Creation of Algebra:";
    AttachSpec("~/packages_github/AlgEt/spec");
    SetAssertions(2);
    _<x>:=PolynomialRing(Integers());
    printf ".";
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    assert #Basis(A) eq Dimension(A);

    seq:=[x^2-5,x^2-7];
    seq:=[NumberField(f) : f in seq];
    A:=EtaleAlgebra(seq);
    printf ".";
    SetAssertions(1);
    printf " all good!\n";

*/
