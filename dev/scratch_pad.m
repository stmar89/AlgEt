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
