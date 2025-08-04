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
// Copyright 2024, S. Marseglia
/////////////////////////////////////////////////////

freeze;

declare attributes AlgEtQ : IsNumberField;

//------------
// Creation for AlgEtQ
//------------

intrinsic IsNumberField(A::AlgEtQ) -> BoolElt,FldNum,Map 
{Given an étale algebra over Q returns wheter it is a number field, and if so the number field and an isomorphism from the étale algebra to the number field.}
   if not assigned A`IsNumberField then
        Ks,embs,projs:=Components(A);
        if #Ks gt 1 then
            return false,_,_;
        else
            K:=Ks[1];
            emb:=embs[1];
            proj:=projs[1];
            isom:=map<A->K| x:->proj(x), y:->emb(y) >;
            return true,K,isom;
        end if;
   end if;
   return Explode(A`IsNumberField);
end intrinsic;

/* TESTS

    printf "### Testing IsNumberField:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    SetAssertions(2);
    _<x>:=PolynomialRing(Integers());
    printf ".";
    f:=x^8+16;
    A:=EtaleAlgebra(f);
    test,K,i:=IsNumberField(A);
    assert test;
    for j in [1..100] do
        a:=Random(A);
        b:=Random(A);
        assert i(a+b) eq i(a)+i(b);
        assert i(a*b) eq i(a)*i(b);
        ii:=Inverse(i);
        OK:=MaximalOrder(K);
        a:=Random(OK,3);
        b:=Random(OK,3);
        assert ii(a+b) eq ii(a)+ii(b);
        assert ii(a*b) eq ii(a)*ii(b);
    end for;



    seq:=[x^2-5,x^2-7];
    seq:=[NumberField(f) : f in seq];
    A:=EtaleAlgebra(seq);
    assert not IsNumberField(A);
    printf ".";

    SetAssertions(1);
    printf " all good!";
*/
