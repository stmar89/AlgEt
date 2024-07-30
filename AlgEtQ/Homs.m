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

declare verbose AlgEtQHoms, 1;

declare attributes AlgEtQ : HomsToC;

import "Ord.m" : crQZ , crZQ , Columns , hnf , MatrixAtoQ , MatrixAtoZ , MatrixQtoA , meet_zbasis , inclusion_matrix;

//------------
// Homomorphism to the Complex Numbers
//------------

intrinsic HomsToC(A::AlgEtQ : Prec:=30)->SeqEnum[Map]
{Returns the sequence of homomorphisms from A to the complex field CC. The precision of CC is given by the optional parameter "Prec". Default value is 30}
    if not assigned A`HomsToC or (Prec ne Precision(Codomain(A`HomsToC[1]))) then
        CC:=ComplexField(Prec);
        images:=function(x)
            return &cat[[CC ! z : z in Conjugates(y : Precision:=Prec)] : y in Components(x)];
        end function;
        maps:=< map< A -> CC | x:-> images(x)[k] > : k in [1..Dimension(A)] >;
        A`HomsToC:=maps;
    end if;
    return A`HomsToC;
end intrinsic;

//------------
// Homomorphism between AlgEtQ
//------------

intrinsic Hom(A::AlgEtQ , B::AlgEtQ , img::SeqEnum[AlgEtQElt] : CheckMultiplicative:=true, CheckUnital:=false, ComputeInverse:=true)->Map
{Given two étale algebras A and B and a sequence img of elements of B, returns the Q-algbra homomorphism defined by sending the AbsoluteBasis of A to img. The VarArg CheckMultiplicative determines if the multiplicativity of the defined map is checked, while the VarArg CheckUnital determines wheter One(A) is sent to One(B). If the VarArg ComputeInverse is true, it checkes whether the map is invertible and, if so, it defines also the inverse (by assigning preimages).}
    basis:=AbsoluteBasis(A);
    require forall{x : x in img | x in B} and #img eq #basis : "the images do not defined an additive map.";
    image:=function(x)
        return SumOfProducts(AbsoluteCoordinates(x),img);
    end function;
    
    has_inverse:=false;
    if ComputeInverse then 
        mm:=MatrixAtoQ(img);
        has_inverse:=NumberOfRows(mm) eq NumberOfColumns(mm) and Determinant(mm) ne 0;
    end if;

    if has_inverse then
        img_inv:=MatrixQtoA(B,mm^-1);
        preimage:=function(y)
            return SumOfProducts(AbsoluteCoordinates(y),img_inv);
        end function;
        m:=map<A->B | x:->image(x), y:->preimage(y)>;
    else
        m:=map<A->B | x:->image(x)>;
    end if;

    if CheckMultiplicative then
        require forall{ x : x,y in basis | m(x)*m(y) eq m(x*y) } : "the images do not define a multiplicative map.";
    end if;

    if CheckUnital then 
        require m(One(A)) eq One(B) : "One(A) is not sent to One(B)";
    end if;
    return m;
end intrinsic;



/* TESTS

    printf "### Testing Homs:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    homs:=HomsToC(A);
    a:=PrimitiveElement(A);
    assert &and[ Abs(Evaluate(f,h(a))) lt 10^-20 : h in homs ];
    old_prec:=Precision(Codomain(homs[1]));
    new_prec:=10*old_prec;
    homs:=HomsToC(A : Prec:=new_prec);
    assert Precision(Codomain(homs[1])) eq new_prec;
    B:=EtaleAlgebra(Components(A) cat [NumberField(x^2+2)]);
    img:=[ B!(Components(b) cat <0>) : b in AbsoluteBasis(A) ];
    incl:=Hom(A,B,img : CheckMultiplicative:=true );
    assert incl(One(A)) ne One(B);
    assert forall{ a : a in AbsoluteBasis(A) | MinimalPolynomial(incl(a)) eq MinimalPolynomial(a)};
    aut:=[ Automorphisms(K) : K in Components(B) ];
    aut:=[ Random(a) : a in aut ];
    img:=[ B!<aut[i](Components(a)[i]) : i in [1..#aut] > : a in AbsoluteBasis(B) ];
    aut:=Hom(B,B,img: CheckUnital:=true);
    inv:=Inverse(aut);
    assert forall{ b : b in AbsoluteBasis(B) | (inv*aut)(b) eq b and (aut*inv)(b) eq b};
    printf " all good!\n"; 

*/
