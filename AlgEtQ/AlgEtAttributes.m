/* vim: set syntax=magma :*/

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, stefano.marseglia89@gmail.com
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

//------------
// Printing
//------------

intrinsic Print(A::AlgEtQ)
{Prints the defining polynomial or the components defining A.}
    if assigned A`DefiningPolynomial then
        f:=DefiningPolynomial(A);
        printf "Etale Algebra over (%o) defined by %o", FieldOfFractions(CoefficientRing(f)),f;    
    else
        printf "Etale Algebra product of %o", Components(A);
    end if;
end intrinsic;

//------------
// Access attributes
//------------

intrinsic DefiningPolynomial(A::AlgEtQ) -> RngUPolElt
{Returns the defining polynomial of A, if the corresponding number fields are distinct.}
    if not assigned A`DefiningPolynomial then
        polys:=[DefiningPolynomial(L) : L in Components(A)];
        require #polys eq #Seqset(polys) : "The number fields defining A are not distinct.";
        A`DefiningPolynomial:=&*(polys);
    end if;
    return A`DefiningPolynomial;
end intrinsic;

intrinsic Components(A::AlgEtQ) -> SeqEnum
{Returns the number fields of which A is a product of,together with embeddings and projections.}
    return Explode(A`Components);
end intrinsic;

intrinsic Dimension(A::AlgEtQ)->RngInt
{Dimension of A.}    
    if not assigned A`Dimension then
        nf:=Components(A);
        require HasBaseField(A) : "The number fields of A shoud all be defined over the same base ring.";
        A`Dimension:=&+[Degree(E) : E in nf];
    end if;
    return A`Dimension;
end intrinsic;

intrinsic AbsoluteDimension(A::AlgEtQ)->RngInt
{Dimension of A over the prime field.}    
    if not assigned A`AbsoluteDimension then
        A`AbsoluteDimension:=&+[AbsoluteDegree(E) : E in Components(A)];
    end if;
    return A`AbsoluteDimension;
end intrinsic;

intrinsic HasBaseField(A::AlgEtQ) -> BoolElt,FldNum
{Returns whether A has common base field. If this is the case it returns it.}
    if not assigned A`HasBaseField then
        nf,embs:=Components(A);
        F:=BaseRing(nf[1]);
        A`HasBaseField:=forall{ E : E in nf[2..#nf] | BaseRing(E) eq F };
        if A`HasBaseField then
            diag:=map< F->A | x:-> A!<nf[i]!x : i in [1..#nf]> >;
            A`BaseField:=<F,diag>;
        end if;
    end if;
    return A`HasBaseField;
end intrinsic;

intrinsic BaseField(A::AlgEtQ) -> FldNum
{Returns the common base field of the Algebra, if it exists.}
    if not assigned A`BaseField then
        require HasBaseField(A) : "The number fields should all be defined over the same Base ring/field.";
        // if HasBaseField is true, then it is assiged
    end if;
    return Explode(A`BaseField);
end intrinsic;

intrinsic PrimeField(A::AlgEtQ) -> FldNum
{Returns the prime field of the Algebra.}
    if not assigned A`PrimeField then
        nf:=Components(A);
        F:=PrimeField(nf[1]); //this should always be Rationals()
        A`PrimeField:=F;
    end if;
    return A`PrimeField;
end intrinsic;

//------------
// Equality
//------------

intrinsic 'eq'(A1::AlgEtQ,A2::AlgEtQ) -> BoolElt
{Equality testing.}
   c1,e1,p1:=Components(A1);
   c2,e2,p2:=Components(A2);
   return <c1,e1,p1> eq <c2,e2,p2>;
end intrinsic;

/* TESTS

    printf "### Testing Attributes and Equality:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    SetAssertions(2);
    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    B:=EtaleAlgebra(f);
    assert not A eq B;
    assert A ne B;
    assert A eq A;
    assert A cmpeq A;
    assert not A cmpeq B;
    SetAssertions(1);
    printf " all good!\n";

*/
