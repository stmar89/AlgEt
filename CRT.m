/* vim: set syntax=magma :*/

//freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose CRT, 2;

declare attributes AlgEtOrd : CRT_data;
declare attributes AlgEtIdl : CRT_data;

/*TODO:

*/

import "Ord.m" : crQZ , crZQ , Columns , hnf , MatrixAtoQ , MatrixAtoZ , MatrixQtoA , meet_zbasis ;

//------------
// Chinese RemainderTheorem
//------------

// // OLD VERSION to be cleaned
// CRT_data_order:=function(S)
//     if not assigned S`CRT_data then
//         K:=Algebra(S);
//         zb:=ZBasis(S); 
//         //I need to modify the ZBasis(S) in a way that One(K) is the first element
//         pos:=Position(zb,One(K));
//         if pos ne 1 then
//             vprintf CRT,2 : "1 is not in the first position of ZBasis(S)\n";
//             if pos eq 0 then //One(K) not in Zbasis_S
//                 vprintf CRT,2 : "1 is not in ZBasis(S)\n";
//                 coord:=AbsoluteCoordinates([One(K)],zb)[1];
//                 // since One(K) is in S, and there exists a Z-Basis of S containing One(K), then there must be a coefficient of coord that is a unit in Z, i.e. 1 or -1.
//                 // the element for which this occurs can be replaced by One(K)
//                 pos:=Position(coord,1);
//                 if pos eq 0 then
//                     pos:=Position(coord,-1);
//                 end if;
//                 assert pos ne 0;
//             end if; 
//             //replacing zb[pos] with One(K) since they generate the same Z-Span
//             temp:=zb[1];
//             zb[1]:=One(K);
//             if pos ne 1 then
//                 zb[pos]:=temp;
//             end if;
//         end if;
//         assert2 Order(zb) eq S;
//         M:=MatrixAtoQ(zb);
//         Minv:=M^-1;
//         S`CRT_data:=<zb,Minv>;
//     end if;
//     return Explode(S`CRT_data);
// end function;
// 
// CRT_data_ideal:=function(I)
//     if not assigned I`CRT_data then
//         _,Minv:=CRT_data_order(Order(I)); 
//         vprintf CRT,2 : "CRT_data_ideal output = %o\n",(MatrixAtoQ(ZBasis(I))*Minv);
//         I`CRT_data:=crQZ(MatrixAtoQ(ZBasis(I))*Minv);
//     end if;
//     return I`CRT_data;
// end function;
// 
// intrinsic ChineseRemainderTheorem(I::AlgEtIdl,J::AlgEtIdl,a::AlgEtElt,b::AlgEtElt)-> AlgEtElt
// {Given two coprime ideals I and J of S, two elements a,b in S, finds e such that (e-a) in I and (e-b) in J.}
//     S:=Order(I);
//     vprintf CRT,2 : "S :=%o;\nI := %o;\nJ := %o;\n\/\/[a,b] =\nelts := %o;\n",PrintSeqAlgEtElt(ZBasis(S)),PrintSeqAlgEtElt(ZBasis(I)),PrintSeqAlgEtElt(ZBasis(J)),PrintSeqAlgEtElt([a,b]); 
//     require a in S and b in S:"the elements must lie in order of definition of the ideals";
//     require S eq Order(J): "the ideals must be of the same order";
//     require IsCoprime(I,J) : "the ideals must be coprime";
//     I_min:=MinimalInteger(I);
//     J_min:=MinimalInteger(J);
//     g,c1,d1:=XGCD(I_min,J_min);
//     if g ne 1 then
//         K:=Algebra(S);
//         n:=AbsoluteDimension(K);
//         Zbasis_S,Minv:=CRT_data_order(S);
//         vprintf CRT,2 : "ZBasis_S := %o;\n",PrintSeqAlgEtElt(Zbasis_S);
//         A:=CRT_data_ideal(I);
//         B:=CRT_data_ideal(J);
//         C:=VerticalJoin(A,B);
//         H,U:=HermiteForm(C); //U*C = H;
//         z:=ZeroMatrix(Integers(),n,n);
//         s:=ScalarMatrix(n,1);
//         assert2 H eq VerticalJoin(s,z);
//         P:=VerticalJoin(HorizontalJoin(z,s),HorizontalJoin(s,z));
//         U1:=Transpose(U)*P; //I need the (n+1)st column of U1
//         Z:=Transpose(U1)[n+1];
//         X:=Matrix(Integers(),1,n,[Z[i] : i in [1..n]]);
//         Y:=X*A;
//         c:=&+[Y[1,i]*Zbasis_S[i] : i in [1..n]];
//         assert2 c in I;
//         d:=One(K)-c;
//         assert2 d in J;
//         // c in I, d in J st 1 = c + d
//     else
//         //g:=c1*I_min+d1*J_min
//         c:=c1*I_min;
//         d:=d1*J_min;
//     end if;
//     e:=a*d+b*c;
//     vprintf CRT,2 : "e := %o;\n",PrintSeqAlgEtElt([e])[1];
//     assert e-a in I;
//     assert e-b in J;
//     return e;
// end intrinsic;

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

intrinsic ChineseRemainderTheorem(Is::SeqEnum[AlgEtIdl],as::SeqEnum[AlgEtElt])-> AlgEtElt
{Given a sequence Is of coprime ideals of S, and a sequence as of elements of S, it returns an element e such that e-as[i] in Is[i] for every i.}
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
        // the following 2 lines are the slowest, many coercions in A
        //cs:=[ &+[cc[i][j]*ZBasis(Is[i])[j] : j in [1..n]] : i in [1..N] ]; 
        cs:=[ SumOfProducts(cc[i],ZBasis(Is[i])) : i in [1..N] ]; 
        // cs[i] in Is[i] and \sum_i cs[i] = 1
    else
        //1 = g = \sum_i c1s[i]*Is_min[i]
        cs:=[c1s[i]*Is_min[i] : i in [1..N]];// only integers here. very fast
        // the following lines is the slowest, many coercions in A
    end if;
    //e:=&+[cs[i]*ashat[i] : i in [1..N]];
    e:=SumOfProducts(cs,ashat);
    assert2 &+cs eq One(Algebra(S));
    assert2 forall{ cs : i in [1..N] | cs[i] in Is[i] };
    vprintf CRT,2 : "e := %o;\n",PrintSeqAlgEtElt([e])[1];
    assert forall{ i : i in [1..N] | e-as[i] in Is[i]};
    return e;
end intrinsic;

intrinsic ChineseRemainderTheorem(I::AlgEtIdl,J::AlgEtIdl,a::AlgEtElt,b::AlgEtElt)-> AlgEtElt
{Given two coprime ideals I and J of S, two elements a,b in S, finds e such that (e-a) in I and (e-b) in J.}
    return ChineseRemainderTheorem([I,J],[a,b]);
end intrinsic;

/* TEST

    AttachSpec("~/packages_github/AlgEt/spec");
    SetVerbose("CRT",1);
    SetAssertions(1);

    ////////////////
    //Simple extension of Q
    ///////////////

    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    E1:=EquationOrder(A);
    
    time pp:=PrimesAbove(Conductor(E1));
    time pp13:=[ P : P in pp | MinimalInteger(P) eq 13 ];

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
    // test 1
    
    t0:=Cputime();
    out1:=[];
    for pair in pairs do
        a:=pair[1];
        b:=pair[2];
        e:=ChineseRemainderTheorem(pp13[1],pp13[2],a,b);
        Append(~out1,e);
    end for;
    Cputime(t0); // old code ~14 secs. new code ~7 secs. w/o profiler

    // test 2
    t0:=Cputime();
    out2:=[];
    e1:=ChineseRemainderTheorem(pp13[1],pp13[2],A!1,A!0);
    e2:=ChineseRemainderTheorem(pp13[1],pp13[2],A!0,A!1);
    for pair in pairs do
        a:=pair[1];
        b:=pair[2];
        e:=a*e1+b*e2;
        Append(~out2,e);
    end for;
    Cputime(t0);
    pp:=pp13[1]*pp13[2];
    assert forall{i : i in [1..#out1] | (out1[i] - out2[i]) in pp};
    
    //////////////////
    //Relative setting
    /////////////////
    // issues here :-(

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


*/
