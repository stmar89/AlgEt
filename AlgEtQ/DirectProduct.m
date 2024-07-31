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
// Copyright 2024, S. Marseglia
/////////////////////////////////////////////////////

freeze;

declare verbose AlgEtQ, 1;

//------------
// DirectProduct of AlgEtQ
//------------

intrinsic DirectProduct(seq::SeqEnum[AlgEtQ]) -> AlgEtQ,SeqEnum[Map],SeqEnum[Map]
{Given a sequence of étale algebras, it returns the direct product,togheter with canonical inclusions and projections.}
    all_comps:=[ Components(K) : K in seq ];
    tot_num_comps:=&+[#c : c in all_comps];
    D:=EtaleAlgebra(&cat(all_comps));
    embs:=[* *];
    projs:=[* *];
    counter:=0;
    for iK->K in seq do
        counter_new:=counter+#all_comps[iK];
        emb_img:=[ D!< counter+1 le j and j le counter_new select Components(a)[j-counter] else 0 : j in [1..tot_num_comps] > : a in AbsoluteBasis(K) ];
        emb:=Hom(K,D,emb_img:CheckMultiplicative:=true, ComputeInverse:=false);
        proj_img:=[ K!< Components(d)[i] : i in [counter+1..counter_new]> : d in AbsoluteBasis(D) ];
        proj:=Hom(D,K,proj_img:CheckMultiplicative:=true, ComputeInverse:=false);
        assert2 forall{ a : a in AbsoluteBasis(K) | proj(emb(a)) eq a };
        counter:=counter_new;
        Append(~embs,emb);
        Append(~projs,proj);
    end for;
    return D,embs,projs;
end intrinsic;

/* TESTS

    printf "### Testing DirectProduct:";
    //AttachSpec("~/AlgEt/spec");
    SetAssertions(2);
    _<x>:=PolynomialRing(Integers());
    seq:=[x^2-5,x^2-7,x^2-5,x^2-11,x^3+x+1];
    seq:=[NumberField(f) : f in seq];
    A:=EtaleAlgebra(seq);
    _,_,_:=DirectProduct([A,A,A]);
    printf ".";

    B:=EtaleAlgebra([NumberField(x^8+16),NumberField(x^3+x+1)]);
    _,_,_:=DirectProduct([A,B]);
    printf ".";
    
    SetAssertions(1);
    printf " all good!\n";

*/
