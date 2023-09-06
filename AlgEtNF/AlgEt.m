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

//freeze;

declare verbose AlgEtNF, 1;

/*TODO:

*/

declare type AlgEtNF[AlgEtNFElt];

declare attributes AlgEtNF : DefiningPolynomial, 
                           // ass_algebra, 
                           Dimension,
                           AbsoluteDimension,
                           BaseField, //a tup : <F,m> where F is the Base field and m is the diagonal embedding into A
                           HasBaseField, //a boolean
                           PrimeField,
                           Components; //a tup of 3 sequences: the first are the NF, 
                                         //the second are embeddings and the third are projections

//------------
// Creation for AlgEtNF
//------------

intrinsic EtaleAlgebra(seq::SeqEnum[FldNum[FldNum]]) -> AlgEtNF
{Given a sequence of number fields returns the étale algebra corresponding to the direct product.}
    A:=New(AlgEtNF);
    embs:=[ map< seq[i]->A | x:-> A! (<seq[j]!0 : j in [1..i-1]> cat <x> cat <seq[j]!0 : j in [i+1..#seq]>)  > : i in [1..#seq] ];
    projs:=[ map< A->seq[i] | y:-> Components(y)[i] > : i in [1..#seq] ];
    A`Components:=<seq,embs,projs>;
    return A;
end intrinsic;

intrinsic EtaleAlgebra(f::RngUPolElt[FldNum]) -> AlgEtNF
{Given a squarefree polynomial over the integers returns the product of the number fields defined by the irreducible factors.}
    require IsSquarefree(f) : "The polynomial must be squarefree.";
    A:=EtaleAlgebra([NumberField(g[1]) : g in Factorization(f)]);
    A`DefiningPolynomial:=f;
    return A;
end intrinsic;

intrinsic EtaleAlgebra(f::RngUPolElt[FldNum]) -> AlgEtNF
{Given a squarefree polynomial over the rationals returns the product of the number fields defined by the irreducible factors.}
    require IsSquarefree(f) : "The polynomial must be squarefree.";
    A:=EtaleAlgebra([NumberField(g[1]) : g in Factorization(f)]);
    A`DefiningPolynomial:=f;
    return A;
end intrinsic;

/* TESTS

    //////////////////
    // Attributes Relative setting
    /////////////////

    // Tests for relative extensions.
    K:=NumberField(x^2-5);
    _<y>:=PolynomialRing(K);
    p:=y^2-7;
    A:=EtaleAlgebra(p);
    F,mFA:=BaseField(A);
    _:=[ mFA(b) : b in Basis(MaximalOrder(F)) ];
    printf ".";

    E:=NumberField(p);
    seq:=[E,K];
    A:=EtaleAlgebra(seq);
    assert not HasBaseField(A); // A is a product of an absolute and a relative extensions
    assert AbsoluteDimension(A) eq 6;
    assert PrimeField(A) eq RationalField();
    printf ".";

    A:=EtaleAlgebra([E,E]);
    assert Dimension(A) eq 4;
    assert AbsoluteDimension(A) eq 8;
    assert BaseField(A) eq K;
    assert PrimeField(A) eq RationalField();
    SetAssertions(1);
    printf ".";
    printf " all good!\n";

    //////////////////
    // CRT Relative setting
    /////////////////

    _<x>:=PolynomialRing(Integers());
    K:=NumberField(x^2-5);
    _<y>:=PolynomialRing(K);
    K1:=NumberField(y^2-49*7*K.1);
    K2:=NumberField(y^5-25*7*K.1);
    A:=EtaleAlgebra([K1,K2]); 
    idem:=OrthogonalIdempotents(A);
    time O:=MaximalOrder(A);
    repeat 
        E:=Order( [Random(O) : i in [1..3]] );
    until exists{ i : i in idem | not i in E} and not 1 in ZBasis(E);
    n:=100;
    repeat
        n+:=1;
        I:=n*E;
    until IsPrime(n) and #PrimesAbove(I) gt 1;
    pp:=PrimesAbove(I);

    pairs:=[];
    for i in [1..100] do
        repeat
            a:=Random(E);
        until not a in pp[1];
        repeat
            b:=Random(E);
        until not b in pp[2];
        Append(~pairs,[a,b]);
    end for;
    t0:=Cputime();
    for pair in pairs do
        e:=ChineseRemainderTheorem(pp[1],pp[2],a,b);
    end for;
    Cputime(t0);


    ///////////////////////
    // Elements Relative Extensions
    ///////////////////////
    K:=NumberField(x^2-5);
    _<y>:=PolynomialRing(K);
    p:=y^2-7;
    A:=EtaleAlgebra(p);
    a:=PrimitiveElement(A);
    assert MinimalPolynomial(a) eq p;
    A!1 - A!(1/2);
    A!1 / A!(1/2);

    
    A:=EtaleAlgebra([K,K]);
    e:=A!<K.1,K.1^2>;
    MinimalPolynomial(e);
    IsIntegral(e);
    e2:=A!<K.1^2,K.1^2>;
    assert MinimalPolynomial(e2) eq MinimalPolynomial(K.1^2);
    _,embs,projs:=Components(A);
    assert embs[1](K.1)+embs[2](K.1^2) eq e;
    assert projs[1](e) eq K.1;
    assert projs[2](e) eq K.1^2;

    E:=NumberField(p);
    seq:=[E,K];
    A:=EtaleAlgebra(seq);
    A!1 + A!(1/2);
    A!<seq[1]!1,seq[2]!(1/2)>/A![1,2];
    assert #AbsoluteBasis(A) eq AbsoluteDimension(A);
    e:=A!<K.1+E.1,K.1^2>;
    assert #AbsoluteCoordinates(e) eq AbsoluteDimension(A);
    assert e eq &+[AbsoluteCoordinates(e)[i]*AbsoluteBasis(A)[i] : i in [1..AbsoluteDimension(A)]];
    OrthogonalIdempotents(A);
    Idempotents(A);

    A:=EtaleAlgebra([E,E]);
    assert #Basis(A) eq Dimension(A);
    assert #AbsoluteBasis(A) eq AbsoluteDimension(A);
    OrthogonalIdempotents(A);
    Idempotents(A);

    ///////////////////////
    // FactPrimes Relative Extensions
    ///////////////////////
    K:=NumberField(x^2-5);
    _<y>:=PolynomialRing(K);
    E1:=NumberField(y^2-49*7*K.1);
    E2:=NumberField(y^5-25*7*K.1);
    A:=EtaleAlgebra([E1,E2]); 
    time O:=MaximalOrder(A);
    time assert IsBass(O);
    time ids:=[ Ideal(O,[Random(O) : i in [1..10]]) : i in [1..100]]; 
    time facs:=[ Factorization(I) : I in ids | I ne OneIdeal(O) ];

    repeat 
        E:=Order( [Random(O) : i in [1..3]] );
    until not IsMaximal(E);
    // _:=IsGorenstein(E); // This triggers an ERROR since TraceDualIdeal is implemented only for AlgEtNF over Q.
    _:=IsBass(E);
    _:=PrimesAbove(Conductor(E));

        ///////////////////////
    // Idl: Relative Extensions
    ///////////////////////
    K:=NumberField(x^2-5);
    _<y>:=PolynomialRing(K);
    E1:=NumberField(y^2-49*7*K.1);
    E2:=NumberField(y^5-25*7*K.1);
    A:=EtaleAlgebra([E1,E2]); 
    assert HasBaseField(A);
    time O:=MaximalOrder(A);
    for i in [1..100] do
        a:=Random(O);
        b:=Random(O);
        I:=Ideal(O,a);
        I:=a*O;
        I:=O*a;
    end for;
    time _:={ Random(O)*O : i in [1..100] };
    _:=a*O + O!!(a*O);
    time _:=&*[Random(O)*O : i in [1..100]];
    time _:=&*[(Random(O)*O+Random(O)*O) : i in [1..100]];
    time I:=Ideal(O,[Random(O) : i in [1..100]]);
    time J:=&*[I : i in [1..100]];
    time JJ:=I^100;
    time assert J eq JJ;
    time _:=&meet[(Random(O)*O+Random(O)*O) : i in [1..100]];
    time assert forall{ I : I in [Random(O)*O : i in [1..100]] | IsProductOfIdeals(I)};


    //////////////////////
    // Ord :  Relative extensions
    //////////////////////
    
    K:=NumberField(x^2-5);
    _<y>:=PolynomialRing(K);
    p:=y^2-7;
    assert forall{z : z in ZBasis(O1) | z in O1 };
    A:=EtaleAlgebra(p);
    OA:=MaximalOrder(A);
    time O1:=Order(Basis(A)); //this should trigger an error
    time O2:=Order(AbsoluteBasis(A));
    assert forall{z : z in ZBasis(O2) | z in O2 };
    for O in [O2] do
        for i in [1..100] do
            assert Random(O) in O;
        end for;
    end for;

    
    
    A:=EtaleAlgebra([K,K]);
    time O1:=Order(Basis(A));
    time O2:=Order(AbsoluteBasis(A));
    time O1 eq O2;
    assert forall{z : z in ZBasis(O1) | z in O1 };
    assert forall{z : z in ZBasis(O2) | z in O2 };
    for O in [O1,O2] do
        for i in [1..100] do
            assert Random(O) in O;
        end for;
    end for;


    seq:=[NumberField(p),NumberField(x^2-5)];
    A:=EtaleAlgebra(seq);
    time O2:=Order(AbsoluteBasis(A));
    for O in [O2] do
        for i in [1..100] do
            assert Random(O) in O;
        end for;
    end for;


    K:=NumberField(x^2-25*5);
    _<y>:=PolynomialRing(K);
    p:=y^2-7;
    E:=NumberField(p); //relative ext
    seq:=[E,E];
    A:=EtaleAlgebra(seq);
    time O2:=Order(AbsoluteBasis(A));
    F,m:=BaseField(A);
    assert F eq K;
    time O1:=Order(Basis(A));//this should trigger an error
    time O3:=OrderOver(Basis(A),EquationOrder(F));
    time O4:=OrderOver(Basis(A),MaximalOrder(F));
    time O3 eq O4;
    time O3 eq O2;
    time O4 eq O2;
    assert forall{ O : O in [O2,O3,O4] | forall{z : z in ZBasis(O) | z in O}};
    for O in [O2,O3,O4] do
        for i in [1..100] do
            assert Random(O) in O;
        end for;
    end for;


    /////////////////////
    // TraceNorm Relative extension
    /////////////////////

    K:=NumberField(x^2-5);
    _<y>:=PolynomialRing(K);
    E1:=NumberField(y^2-49*7*K.1);
    E2:=NumberField(y^5-25*7*K.1);
    A:=EtaleAlgebra([E1,E2]); 
    assert HasBaseField(A);
    O:=MaximalOrder(A);
    for i in [1..100] do
        a:=Random(O);
        b:=Random(O);
        assert AbsoluteTrace(a)+AbsoluteTrace(b) eq AbsoluteTrace(a+b);
        assert Trace(a)+Trace(b) eq Trace(a+b);
        assert AbsoluteNorm(a)*AbsoluteNorm(b) eq AbsoluteNorm(a*b);
        assert Norm(a)*Norm(b) eq Norm(a*b);
    end for;
    
    A:=EtaleAlgebra([K,K]);
    for i in [1..100] do
        a:=Random(A);
        b:=Random(A);
        assert Trace(a)+Trace(b) eq Trace(a+b);
        assert Norm(a)*Norm(b) eq Norm(a*b);
    end for;

    K:=NumberField(x^2-25*5);
    _<y>:=PolynomialRing(K);
    p:=y^2-7;
    E:=NumberField(p); //relative ext
    A:=EtaleAlgebra([E,E]);
    O:=MaximalOrder(A);
    for i in [1..100] do
        a:=Random(O);
        b:=Random(O);
        assert AbsoluteTrace(a)+AbsoluteTrace(b) eq AbsoluteTrace(a+b);
        assert Trace(a)+Trace(b) eq Trace(a+b);
        assert AbsoluteNorm(a)*AbsoluteNorm(b) eq AbsoluteNorm(a*b);
        assert Norm(a)*Norm(b) eq Norm(a*b);
    end for;

*/
