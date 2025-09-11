/////////////////////////////////////////////////////
// Copyright 2025.
// Stefano Marseglia, stefano.marseglia89@gmail.com
// https://stmar89.github.io/index.html
//
// Distributed under the terms of the CC-BY 4.0 licence.
// https://creativecommons.org/licenses/by/4.0/
/////////////////////////////////////////////////////



freeze;

declare verbose AlgEtQ, 1;

//------------
// DirectProduct of AlgEtQ
//------------

///## Direct products

/// Given a sequence of étale algebras over $\mathbb{Q}$, returns the direct product, together with sequences of inclusions and projections.
intrinsic DirectProduct(seq::SeqEnum[AlgEtQ]) -> AlgEtQ,SeqEnum[Map],SeqEnum[Map]
{Given a sequence of étale algebras over the rational field, returns the direct product, together with sequences of inclusions and projections.}
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
    printf " all good!";

*/