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

declare verbose AlgEtQ, 1;

declare type AlgEtQ[AlgEtQElt];

declare attributes AlgEtQ : DefiningPolynomial, 
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

intrinsic EtaleAlgebra(seq::SeqEnum[FldNum]) -> AlgEtQ
{Given a sequence of number fields returns the étale algebra corresponding to the direct product. Note: the number fields with DefiningPolynomial of degree one should be created with the vararg DoLinearExtention set to true.}
// funny stuff: ExtendedType([NumberField(g[1]) : g in Factorization(f)]); for f:=(x^8+16)*(x^8+81); returns SeqEnum[FldNum] and not SeqEnum[FldNum[FldRat]]
    require forall{ K : K in seq | ISA(ExtendedType(K),FldNum[FldRat])} : "The given number fields are not absolute extensions of Q.";
    A:=New(AlgEtQ);
    embs:=[ map< seq[i]->A | x:-> A! (<seq[j]!0 : j in [1..i-1]> cat <x> cat <seq[j]!0 : j in [i+1..#seq]>)  > : i in [1..#seq] ];
    projs:=[ map< A->seq[i] | y:-> Components(y)[i] > : i in [1..#seq] ];
    A`Components:=<seq,embs,projs>;
    return A;
end intrinsic;

intrinsic EtaleAlgebra(f::RngUPolElt[RngInt]) -> AlgEtQ
{Given a squarefree polynomial over the integers returns the product of the number fields defined by the irreducible factors.}
    require IsSquarefree(f) : "The polynomial must be squarefree.";
    A:=EtaleAlgebra([NumberField(g[1]: DoLinearExtension) : g in Factorization(f)]); 
    // DoLinearExtension is to avoid issues with g's of degree 1.
    // This can occur while computing TotallyRealSubAlgebra's, for example
    A`DefiningPolynomial:=f;
    return A;
end intrinsic;

intrinsic EtaleAlgebra(f::RngUPolElt[FldRat]) -> AlgEtQ
{Given a squarefree polynomial over the rationals returns the product of the number fields defined by the irreducible factors.}
    require IsSquarefree(f) : "The polynomial must be squarefree.";
    A:=EtaleAlgebra([NumberField(g[1] : DoLinearExtension) : g in Factorization(f)]);
    // DoLinearExtension is to avoid issues with g's of degree 1.
    // This can occur while computing TotallyRealSubAlgebra's, for example
    A`DefiningPolynomial:=f;
    return A;
end intrinsic;

/* TESTS

    printf "### Testing Creation of Algebra:";
    //AttachSpec("~/packages_github/AlgEt/spec");
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

    seq:=[x-1,x-20];
    A:=EtaleAlgebra(&*seq);
    printf ".";

    SetAssertions(1);
    printf " all good!\n";

*/
