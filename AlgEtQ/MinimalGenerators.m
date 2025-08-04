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

//------------
// TwoGeneratingSet
//------------
import "Ord.m" : MatrixAtoQ,MatrixAtoZ;

intrinsic TwoGeneratingSet(I::AlgEtQIdl)
{A procedure that given an invertible ideal I put in the attibute Generators of I two non-zerodivisors in I that generate I. If I is known to be principal nothing is done.}
    require IsInvertible(I) : "the ideal must be invertible";
    if #Generators(I) gt 2 then
        // if assigned I`IsIntegral and IsIntegral(I) then
        //     a:=MinimalInteger(I);
        // else
        //     a:=Random(I : ZeroDivisorsAllowed:=false );
        // end if;
        a:=ShortElement(I);
        S:=Order(I);
        Q,q:=Quotient(I,[a*z : z in ZBasis(S)]);
        if IsTrivial(Q) then
            I`Generators:=[a]; //the ideal is principal
        else
            repeat
                repeat
                    b:=Random(Q);
                until b ne Zero(Q);
                b:=b@@q;
            //until I eq Ideal(S,[a,b]);
            until Q eq sub<Q|[q(b*z):z in ZBasis(S)]>;
            I`Generators:=[a,b];
        end if;
    end if;
end intrinsic;

/* TESTS

    printf "### Testing MinimalGenerators:";
	//AttachSpec("~/packages_github/AlgEt/spec");
    SetClassGroupBounds("GRH");
	_<x>:=PolynomialRing(Integers());
    f:=x^4-1000*x^3-1000*x^2-1000*x-1000;
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    P,p:=PicardGroup(E : GRH:=true); //~10 secs

    for g in Generators(P) do 
        I:=p(g);
        TwoGeneratingSet(I);
        assert #Generators(I) le 2;
        printf ".";
    end for;


    // test if TwoGeneratingSet makes the power faster
    // Conlcusion: yes. By quite a bit!
    f:=x^4-100*x^3-100*x^2-100*x-100;
    A:=EtaleAlgebra(f);
	E:=EquationOrder(A);
    P,p:=PicardGroup(E : GRH:=true);
    repeat
        Ii:=Random(P);
    until Ii ne Zero(P);
    I:=p(Ii);

    delete I`IsInvertible;
    exp:=[ Random(2,30) : i in [1..100]];
    l1:=[ I^i : i in exp ];
    printf ".";

    assert IsInvertible(I);
    TwoGeneratingSet(I);
    assert #Generators(I) eq 2;
    l2:=[ I^i : i in exp ];
    assert l1 eq l2;
    printf ".";

    I:=SmallRepresentative(I);
    delete I`IsInvertible;
    l1:=[ I^i : i in exp ];

    assert IsInvertible(I);
    TwoGeneratingSet(I);
    l2:=[ I^i : i in exp ];
    assert l1 eq l2;
    printf " all good!\n"; 

*/
