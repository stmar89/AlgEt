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

declare verbose CRT, 2;

declare attributes AlgEtQOrd : CRT_data;
declare attributes AlgEtQIdl : CRT_data;

import "Ord.m" : crQZ , crZQ , Columns , hnf , MatrixAtoQ , MatrixAtoZ , MatrixQtoA , meet_zbasis ;

CRT_data_order:=function(S)
    if not assigned S`CRT_data then
        M:=MatrixAtoQ(ZBasis(S));
        d:=Denominator(M);
        _,V:=HermiteForm(crQZ(d*M)); // V*d*M eq H
        cc:=Matrix([AbsoluteCoordinates([One(Algebra(S))],S)[1]])*crZQ(V^-1);
        S`CRT_data:=<cc,d>;
    end if;
    return Explode(S`CRT_data);
end function;

CRT_data_ideal:=function(I)
     if not assigned I`CRT_data then
         _,d:=CRT_data_order(Order(I));
         I`CRT_data:=crQZ(d*MatrixAtoQ(ZBasis(I)));
     end if;
     return I`CRT_data;
end function;

intrinsic ChineseRemainderTheorem(I::AlgEtQIdl,J::AlgEtQIdl,a::AlgEtQElt,b::AlgEtQElt)-> AlgEtQElt
{Given two coprime ideals I and J of S, two elements a,b in S, finds e such that (e-a) in I and (e-b) in J.}
    Is:=[I,J];
    as:=[a,b];
    // the following code works only if N =2
    N:=#as;
    S:=Order(Is[1]);
    require #Is eq N: "The number of ideals is not the same as the number of elements";
    require forall{i : i in [1..N] | as[i] in S}:"the elements must lie in order of definition of the ideals";
    require forall{i : i in [2..N] | Order(Is[i]) eq S}:"the ideals must be of the same order";
    ashat:=[ &*[as[j] : j in [1..N] | i ne j ]  : i in [1..N] ];
    Is_min:=[ MinimalInteger(I) : I in Is ];
    g,c1s:=XGCD(Is_min);
    if g ne 1 then
        // cc=coord of 1 in S * V^-1, where V*(zb(S)) = H. H is the same as below.
        // let U*VerticalJoin(zb(I1),...,zb(IN)) = H 
        // extend cc with zeros.
        // hence cc*U gives the coordinates of 1 in I1+I2+...+IN. 
        // Note that we use only the first n rows of U.
        // So instead of extending cc, we can replace U wit only its first n rows.
        K:=Algebra(S);
        n:=AbsoluteDimension(K);
        cc,d:=CRT_data_order(S);
        C:=VerticalJoin([ CRT_data_ideal(I) : I in Is ]);
        H,U:=HermiteForm(C); //U*C = H;
        cc:=cc*crZQ(Matrix(Rows(U)[1..n]));
        cc:=Partition(Eltseq(cc),n);
        //cs:=[ &+[cc[i][j]*ZBasis(Is[i])[j] : j in [1..n]] : i in [1..N] ]; 
        cs:=[ SumOfProducts(cc[i],ZBasis(Is[i])) : i in [1..N] ]; 
        // cs[i] in Is[i] and \sum_i cs[i] = 1
    else
        //1 = g = \sum_i c1s[i]*Is_min[i]
        cs:=[c1s[i]*Is_min[i] : i in [1..N]];// only integers here. very fast
    end if;
    //e:=&+[cs[i]*ashat[i] : i in [1..N]];
    e:=SumOfProducts(cs,ashat);
    assert2 &+cs eq One(Algebra(S));
    assert2 forall{ cs : i in [1..N] | cs[i] in Is[i] };
    assert forall{ i : i in [1..N] | e-as[i] in Is[i]};
    return e;
end intrinsic;

intrinsic ChineseRemainderTheorem(Is::SeqEnum[AlgEtQIdl],as::SeqEnum[AlgEtQElt])-> AlgEtQElt
{Given a sequence `Is` of ideals of S, pairwise coprime, and a sequence `as` of elements of S, it returns an element e such that e-as[i] in Is[i] for every i.}
    N:=#as;
    S:=Order(Is[1]);
    require #Is eq N: "The number of ideals is not the same as the number of elements";
    require forall{i : i in [1..N] | as[i] in S}:"the elements must lie in order of definition of the ideals";
    require forall{i : i in [2..N] | Order(Is[i]) eq S}:"the ideals must be of the same order";
    if N eq 1 then 
        return as[1];
    else
        b:=$$(Is[2..N],as[2..N]);
        return ChineseRemainderTheorem(Is[1],&meet(Is[2..N]),as[1],b);
    end if;
end intrinsic;

intrinsic ChineseRemainderTheoremFunctions(Is::SeqEnum[AlgEtQIdl])-> AlgEtQElt
{Given a sequence `Is` of N ideals of S, pairwise coprime, returns a function S->S^N representing the natural isomorphism S/&*(Is) -> \prod_(I in Is) S/I and a function S^N-S representing the inverse.}
    S:=Order(Is[1]);
    N:=#Is;
    require forall{i : i in [2..N] | Order(Is[i]) eq S}:"the ideals must be of the same order";
    S:=OneIdeal(S);
    Q,q:=Quotient(S,&meet(Is));
    quots:=[];
    maps:=<>;
    for I in Is do
        QI,qI:=Quotient(S,I);
        Append(~quots,QI); 
        Append(~maps,qI); 
    end for;
    D,embs,projs:=DirectSum(quots);
    assert IsIsomorphic(D,Q);
    isom:=iso<Q->D | [ &+[embs[j](maps[j](Q.i@@q)) : j in [1..#Is]] : i in [1..Ngens(Q)] ]>;
    // isom : E1/&meet(Is) -> \prod_{I in Is} E1/I
    // is the natural isomorphism of E1 modules. 
    // The inverse (constructed by considering isom as an addive map) is automatically E1 linear
    func1:=function(x)
        return [projs[j](isom(q(x)))@@maps[j] : j in [1..N] ];
    end function;
    func2:=function(as)
        assert #as eq N;
        return (&+[embs[j](maps[j](as[j])) : j in [1..N] ])@@isom@@q;
    end function;
    II:=&meet(Is);
    assert forall{s : s in ZBasis(S) | func2(func1(s)) -s in II};
    return func1,func2;
end intrinsic;


/* TESTS

    printf "### Testing CRT:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    SetVerbose("CRT",1);
    SetAssertions(2);

    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    E1:=EquationOrder(A);
    
    pp:=PrimesAbove(Conductor(E1));
    pp13:=[ P : P in pp | MinimalInteger(P) eq 13 ];

    pairs:=[];
    for i in [1..10000] do
        repeat
            a:=Random(E1);
        until not a in pp13[1];
        repeat
            b:=Random(E1);
        until not b in pp13[2];
        Append(~pairs,[a,b]);
    end for;
    printf ".";
    // test 1
    
    out1:=[];
    for pair in pairs do
        a:=pair[1];
        b:=pair[2];
        e:=ChineseRemainderTheorem(pp13[1],pp13[2],a,b);
        Append(~out1,e);
    end for;
    printf ".";
    // old code ~14 secs. new code ~7 secs. w/o profiler

    // test 2
    out2:=[];
    e1:=ChineseRemainderTheorem(pp13[1],pp13[2],A!1,A!0);
    e2:=ChineseRemainderTheorem(pp13[1],pp13[2],A!0,A!1);
    for pair in pairs do
        a:=pair[1];
        b:=pair[2];
        e:=a*e1+b*e2;
        Append(~out2,e);
    end for;
    printf ".";
    pp13:=pp13[1]*pp13[2];
    assert forall{i : i in [1..#out1] | (out1[i] - out2[i]) in pp13};

    // test 3 : >2 primes
    tuples:=[]; 
    for i in [1..1000] do
        i_tup:=[];
        for j in [1..#pp] do
            repeat
                a:=Random(E1);
            until not a in pp[j];
            Append(~i_tup,a);
        end for;
        Append(~tuples,i_tup);
    end for;

    out3:=[];
    out4:=[];
    for tup in tuples do
        e:=ChineseRemainderTheorem(pp,tup);
        Append(~out3,e);
    end for;
    
    f,g:=ChineseRemainderTheoremFunctions(pp);
    for tup in tuples do
        e:=g(tup);
        Append(~out4,e);
    end for;

    I:=&meet(pp);
    assert forall{i : i in [1..#out3] | (out3[i] - out4[i]) in I};

    printf ".";

    SetAssertions(1);    
    printf " all good!\n"; 

*/
