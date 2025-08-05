/* vim: set syntax=magma :*/

freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, stefano.marseglia89@gmail.com
// https://stmar89.github.io/index.html
/////////////////////////////////////////////////////

declare verbose ?????, 1;

/*TODO:

*/


//------------
// a so-far-failed attempt to QuotientVS
//------------

// intrinsic QuotientVS(I::AlgEtQIdl, J::AlgEtQIdl, P::AlgEtQIdl) -> ModRng, Map
// {BUGGED Let I, J, P be fractional R-ideals such that:
//  - P is prime of of some order R, with residue field K;
//  - J in I and I/J is a vector space V over K, say of dimension d.
//  The function returns the KModule K^d=V and the natural surjection I->V (with pre-image)}
//     require J subset I : "Teh second argument should be a subset of the first.";
// 	S := Order(P);
//     assert2 P*(Ideal(S,ZBasis(I))) subset Ideal(S,ZBasis(J));
//     assert2 S subset MultiplicatorRing(I);
//     assert2 S subset MultiplicatorRing(J);
//     K,k:=ResidueField(P);
// 	A := Algebra(S);
//     Q,q:=Quotient(I,J);
// 	d := Ilog(#K,#Q); // d = dim(I/J) over (S/P)
//     // SOMETHING WRONG IN HERE. I Think that there are issues when K is not a prime field
//     matrices:=[Matrix(K,[ Eltseq(Q!(q(zS*(Q.j@@q)))) : j in [1..Ngens(Q)]]) : zS in ZBasis(S) ];
//     matrices:=Setseq(Seqset(matrices)); //possibly there are repetiions.
//     V:=RModule(matrices);
//     assert Dimension(V) eq d;
//     assert #Basis(V) eq Ngens(Q);
//     bijVQ:=function(y)
//         eltseq:=Eltseq(y);
//         N:=#Eltseq(eltseq[1]);
//         // y is represented in V wrt to generators that are in bijection with Q. 
//         // So even if a priori the coeficients of each coordinates of y are in the finite field K, 
//         // they are actually all in the prime field of K, and hence can be coerced into integers.
//         assert forall{i : i in [1..Ngens(Q)] | IsZero(Eltseq(eltseq[i])[2..N]) };
//         out:=&+[Q.i*(Integers()!Eltseq(eltseq[i])[1]) : i in [1..Ngens(Q)]];
//         return out; 
//     end function;
//     bij:=map< Q->V | x :-> &+[V.i*Eltseq(x)[i] : i in [1..Ngens(Q)]],
//                      y :-> bijVQ(y) >;
//     v:=map< A-> V | x:->bij(q(x)) , y:->(y@@bij)@@q >;
//     return V,v;
// end intrinsic;

// A new attempt: I was expecting to be faster, but it is not.
intrinsic MaximalIntermediateIdeals_NEW(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]
{Given fractional S-ideals J subseteq I, returns the maximal (with respect to inclusion) fractional S-ideals K such that J subseteq K subsetneq I. Note I is never in the output, while J is in the output if and only if the S-module I/J is simple, in which case the output consists only of J.}
    assert2 J subset I; // "the ideal J needs to be inside I";
    S:=Order(I);
    assert2 S eq Order(J); // "The ideals must be over the same order";
    if J eq I then 
        return {@ @};
    else
        // If J c K c I is maximal then I/K is a simple module, that is, isomorphic to S/P for some prime P of S
        // lying in the support of the finite module I/J.
        // These are precisely the prime P satisfying (J:I) meet S subseteq P.
        //zJ:=ZBasis(J);
        supp:=PrimesAbove(ColonIdeal(J,I) meet OneIdeal(S));
        assert #supp ge 1;
        max_ideals:={@ @};
        for P in supp do
            PIJ:=P*I+J;
            zz:=ZBasis(PIJ);
            Q,q:=QuotientVS(I,PIJ,P);
            max_sub_vs:=MaximalSubmodules(Q);
            assert #max_sub_vs gt 0;
            for V in max_sub_vs do
                K:=Ideal(S,[b@@q:b in Basis(V)] cat zz);
                //assert Index(I,K) eq Index(S,P);
                Include(~max_ideals,K);
            end for;
        end for;
        return max_ideals;
    end if;
end intrinsic;

/////////////
// test &* if the generators are bad it is much better.
/////////////
    //Conlcusion: classic is the winner
    
	AttachSpec("~/packages_github/AlgEt/spec");
	_<x>:=PolynomialRing(Integers());
    f:=x^4-100*x^3-100*x^2-100*x-100;
    //f:=x^4-1000*x^3-1000*x^2-1000*x-1000;
    A:=EtaleAlgebra(f);
	E:=EquationOrder(A);
    time P,p:=PicardGroup(E);
    l:=[];
    for i in [1..10] do
        Ii:=p(Random(P));
        I:=Ii^Random(2,30);
        Append(~l,I);
    end for;

    seq:=l;

    //classic + small rep
    t0:=Cputime();
    seq1:=[];
    as:=[];
    for I in seq do
        aI,a:=SmallRepresentative(I);
        Append(~seq1,aI);
        Append(~as,a);
    end for;
    I3_small:=&*seq1;
    I3:=I3_small*(1/&*as);
    Cputime(t0);

    //classic
    time I4:=&*seq;

    // one creation + small rep
    t0:=Cputime();
    seq1:=[];
    as:=[];
    for I in seq do
        aI,a:=SmallRepresentative(I);
        Append(~seq1,aI);
        Append(~as,a);
    end for;
    gens:=[ Generators(I) : I in seq1 ];
    cc:=CartesianProduct(gens);
    I2_small:=Ideal(E,[&*[d :d in c] : c in cc ]);
    I2:=I2_small*(1/&*as);
    Cputime(t0);

    // one creation : it seems super slow
    t0:=Cputime();
    gens:=[ Generators(I) : I in seq ];
    cc:=CartesianProduct(gens);
    I1:=Ideal(E,[&*[d :d in c] : c in cc ]);
    Cputime(t0);

    assert 1 eq #{I1,I2,I3,I4};
