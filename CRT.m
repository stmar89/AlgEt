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

import "Ord.m" : crQZ , crZQ , Columns , hnf , MatrixAtoQ , MatrixAtoZ , MatrixQToA , meet_zbasis ;

//------------
// Chinese RemainderTheorem
//------------

CRT_data_order:=function(S)
    if not assigned S`CRT_data then
        K:=Algebra(S);
        Zbasis_S:=ZBasis(S);
        //I need to modify the ZBasis(S) in a way that One(K) is the first element of Zbasis_S
        // it is a bit weird, but the output is asserted...
        pos:=Position(Zbasis_S,One(K));
        if pos ne 1 then
            vprintf CRT,2 : "1 is not in the first position of ZBasis(S)\n";
            if pos eq 0 then //One(K) not in Zbasis_S
                // this should happen only if one of the orthogonal idempotens is there.
                // there is no harm at replacing any orthogonal idempotent with One(K).
                vprintf CRT,2 : "1 is not in ZBasis(S)\n";
                coord:=AbsoluteCoordinates([One(K)],Zbasis_S)[1];
                pos:=Position(coord,1);
                if pos eq 0 then
                    pos:=Position(coord,-1);
                end if;
                assert pos ne 0;
            end if; 
            //replacing Zbasis_S[pos] with One(K) since they generate the same Z-Span
            temp:=Zbasis_S[1];
            Zbasis_S[1]:=One(K);
            Zbasis_S[pos]:=temp;
        end if;
        vprintf CRT,2 : "Not using MinimalInteger\n";
        M:=MatrixAtoQ(Zbasis_S);
        // I was expecting this Matrix inversion to be the most expensive one. But apparently it is not.
        Minv:=M^-1;
        S`CRT_data:=<Zbasis_S,Minv>;
    end if;
    return Explode(S`CRT_data);
end function;

CRT_data_ideal:=function(I)
    if not assigned I`CRT_data then
        _,Minv:=CRT_data_order(Order(I)); 
        I`CRT_data:=crQZ(MatrixAtoQ(ZBasis(I))*Minv);
    end if;
    return I`CRT_data;
end function;

intrinsic ChineseRemainderTheorem(I::AlgEtIdl,J::AlgEtIdl,a::AlgEtElt,b::AlgEtElt)-> AlgEtElt
{Given two coprime ideals I and J of S, two elements a,b in S, finds e such that (e-a) in I and (e-b) in J.}
    require IsCoprime(I,J) : "the ideals must be coprime";
    S:=Order(I);
    require a in S and b in S:"the elements must lie in order of definition of the ideals";
    require S eq Order(J): "the ideals must be of the same order";
    I_min:=MinimalInteger(I);
    J_min:=MinimalInteger(J);
    g,c1,d1:=XGCD(I_min,J_min);
    if g ne 1 then
        K:=Algebra(S);
        n:=AbsoluteDimension(K);
        Zbasis_S,Minv:=CRT_data_order(S);
        A:=CRT_data_ideal(I);
        B:=CRT_data_ideal(J);
        C:=VerticalJoin(A,B);
        H,U:=HermiteForm(C); //U*C = H;
        z:=ZeroMatrix(Integers(),n,n);
        s:=ScalarMatrix(n,1);
        assert2 H eq VerticalJoin(s,z);
        P:=VerticalJoin(HorizontalJoin(z,s),HorizontalJoin(s,z));
        U1:=Transpose(U)*P; //I need the (n+1)st column of U1
        Z:=Transpose(U1)[n+1];
        X:=Matrix(Integers(),1,n,[Z[i] : i in [1..n]]);
        Y:=X*A;
        c:=&+[Y[1,i]*Zbasis_S[i] : i in [1..n]];
        assert2 c in I;
        d:=One(K)-c;
        assert2 d in J;
        // c in I, d in J st 1 = c + d
    else
        //g:=c1*I_min+d1*J_min
        c:=c1*I_min;
        d:=d1*J_min;
    end if;
    e:=a*d+b*c;
    assert e-a in I;
    assert e-b in J;
    return e;
end intrinsic;


/* TEST

    Attach("~/packages_github/AlgEt/AlgEt.m");
    Attach("~/packages_github/AlgEt/Elt.m");
    Attach("~/packages_github/AlgEt/Ord.m");
    Attach("~/packages_github/AlgEt/TraceNorm.m");
    Attach("~/packages_github/AlgEt/Idl.m");
    Attach("~/packages_github/AlgEt/WkTesting.m");
    Attach("~/packages_github/AlgEt/FactPrimes.m");
    Attach("~/packages_github/AlgEt/CRT.m");

    SetVerbose("CRT",1);
    SetAssertions(1);

    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    E1:=EquationOrder(A);
    
    time pp:=PrimesAbove(Conductor(E1));
    time pp13:=[ P : P in pp | MinimalInteger(P) eq 13 ];

    pairs:=[];
    for i in [1..100] do
        repeat
            a:=Random(E1);
        until not a in pp13[1];
        repeat
            b:=Random(E1);
        until not b in pp13[2];
        Append(~pairs,[a,b]);
    end for;
    t0:=Cputime();
    SetProfile(true);
    for pair in pairs do
        e:=ChineseRemainderTheorem(pp13[1],pp13[2],a,b);
    end for;
    SetProfile(false);
    ProfilePrintByTotalTime(ProfileGraph());
    Cputime(t0);


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
