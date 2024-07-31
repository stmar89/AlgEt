/* vim: set syntax=magma :*/

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

declare verbose ShortEltSmallRep, 2;
import "Ord.m" : MatrixQtoA,MatrixAtoQ,MatrixAtoZ;

declare attributes AlgEtQIdl : ShortElement, SmallRepresentative;

//------------
// ShortElement
//------------

intrinsic ShortElement(I::AlgEtQIdl) ->AlgEtQElt
{Given an ideal I returns a non-zerodivisor in I with small coefficients (in the LLL sense). This is achieved by randomly picking an element with small coefficients in a LLL-reduced basis (wrt the T2 norm as a Z-lattice).}
    if not assigned I`ShortElement then
        ZBasisLLL(I);
        B := ZBasis(I);
        bd:=1;
        elts:=[ x : x in B | not IsZeroDivisor(x) ];
        while #elts eq 0 do //if all ZeroDiv
            rndm_coeffs:=[ [Random([-bd..bd]) : i in [1..#B]] : j in [1..10]];
            elts:=[ SumOfProducts(rndm,B) : rndm in rndm_coeffs];
            elts:=[ x : x in elts | not IsZeroDivisor(x) ];
            bd +:=1;
        end while;
        _,i:=Min([Norm(x) : x in elts]);
        elt:=elts[i];
        
        // // The following is deterministic, but does not always work. Also, enumerating all short vectors is too much, and too memory extensive.
        // L:=Lattice(MatrixAtoQ(ZBasis(I)));
        // // b:=Basis(LLL(L)); //we reduce above, so it is cached.
        // b:=Basis(L);
        // b:=[ Norm(c) : c in b ];
        // k:=0;
        // stop:=false;
        // min:=Min(b);
        // repeat
        //     k,Round(2^(-k)*min),Round((min+1)*2^k);
        //     p:=ShortVectors(L,2^(-k)*min,(min+1)*2^k);
        //     for i in [1..#p] do
        //         elt:=MatrixQtoA(Algebra(I),Matrix([Eltseq(p[i][1])]))[1];
        //         if not IsZeroDivisor(elt) then
        //             stop:=true;
        //             break i;
        //         end if;
        //     end for;
        //     k+:=1;
        // until stop;
        I`ShortElement:= elt;
    end if;
    return I`ShortElement;
end intrinsic;

//------------
// SmallRepresentative
//------------

intrinsic SmallRepresentative(I::AlgEtQIdl) ->AlgEtQIdl,AlgEtQElt
{Given a fractional R-ideal I, it returns an isomorphic ideal a*I, and the element a, such that a*I is a subset of R, and the cardinality of R/aI is small. This is achieved by computing the ShortElement a of (R:I). Note that if I is invertible R/aI is isomorphic to (R:I)/aR.}
    if not assigned I`SmallRepresentative then                                            
        R:=Order(I);
        cRI:=ColonIdeal(R,I);
        a:=ShortElement(cRI);
        aI:=a*I;
        // the ZBasis of aI might be very big. We make it smaller.
        ZBasisLLL(aI);
        vprintf ShortEltSmallRep,2: "SmallRepresentative:\n
                                I = %o\n,aI = %o\n",PrintSeqAlgEtQElt(ZBasis(I)),PrintSeqAlgEtQElt(ZBasis(aI));
        I`SmallRepresentative:=<aI,a>;
    end if;
    return Explode(I`SmallRepresentative);
end intrinsic;


/* TESTS
    
    printf "### Testing ShortEltSmallRep:";
	//AttachSpec("~/packages_github/AlgEt/spec");
	_<x>:=PolynomialRing(Integers());
    f:=(x^2+5)*(x^2+7)*(x^2+11);
    assert IsSquarefree(f);
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    ff:=Conductor(E);
    _:=ShortElement(ff);
    oo:=FindOverOrders(E); 
    for S in oo do
        printf ".";
        ff:=Conductor(S);
        _:=ShortElement(ff);
    end for;
    printf " all good!\n"; 

*/
