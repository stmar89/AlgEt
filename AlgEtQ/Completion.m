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

intrinsic Completion(P::AlgEtQIdl : MinPrecision:=20) -> FldPad,Map
{Given a prime ideal of the maximal order of an etale algebra L it returns the p-adic field corresponding to the completion LP and a homormophism map:L->LP. The vararg MinPrecision is passed to Completion. map has preimage.}
    L:=Algebra(P);
    require IsMaximal(Order(P)) and IsPrime(P) : "the ideal must be a prime ideal of the maximal order";
    nfs,embs:=Components(L);
    test,PKs:=IsProductOfIdeals(P);
    assert test;
    ind:=[ i : i in [1..#PKs] | not (Order(PKs[i]) ! 1) in PKs[i] ];
    assert #ind eq 1;
    ind:=ind[1];
    K:=nfs[ind];
    mK:=embs[ind]; // mK:=K->L
    PK:=PKs[ind];
    LP,mLP:=Completion(K,PK : Precision:=MinPrecision); // mLP:K->LP
    // it seems that Magma Ignore the Parameter Precision and returns always a ring with Infinity precition.
    // workaround
    LP1:=ChangePrecision(LP,MinPrecision);
    map:=hom< L->LP1 | x:->LP1!mLP(Components(x)[ind]) ,
                      y:->mK((LP!y)@@mLP) >;
    return LP1,map;
end intrinsic;

/* TESTS
    
    printf "### Testing Completion:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    PP<x>:=PolynomialRing(Integers());
    polys:=[
        x^6+3*x^4-10*x^3+15*x^2+125,
        (x^2+5)*(x^4-4*x^3+5*x^2-20*x+25),
        (x^4-5*x^3+15*x^2-25*x+25)*(x^4+5*x^3+15*x^2+25*x+25)
        ];
    for h in polys do
        L:=EtaleAlgebra(h);
        a:=PrimitiveElement(L);
        O:=MaximalOrder(L);
        p:=5;
        pp:=PrimesAbove(p*O);
        for P in pp do
            C,mC:=Completion(P);
        end for;
        printf ".";
    end for;
    printf " all good!\n";

*/
