/* vim: set syntax=magma :*/

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
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
    printf " all good!\n"; 

*/
