/* vim: set syntax=magma :*/

freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose ?????, 1;

/*TODO:

*/


//------------
// What do we do here?
//------------

// NEW idea. it is slower.
quo_frac_idl:=function(I,J)
// I->I/J
    N:=AbsoluteDimension(Algebra(I));
    F:=FreeAbelianGroup(N);
    zbI:=ZBasis(I);
    matI:=MatrixAtoQ(zbI);
    matJ:=MatrixAtoZ(ZBasis(J));
    rel:=[F ! Eltseq(x) : x in Rows(matJ*matI^-1)];
    Q,q:=quo<F|rel>;
    qq:=hom<I -> Q | x:->q(F!AbsoluteCoordinates([x],zbI)[1]) , y:->&+[Eltseq(y@@q)[i]*zbI[i] : i in [1..N]] >;
    return Q,qq;
end function;

intrinsic ChineseRemainderTheorem2(I::AlgEtIdl,J::AlgEtIdl)-> Map,Map,Map
{Given two coprime ideals I and J of S, two elements a,b in S, finds e such that (e-a) in I and (e-b) in J.}
    S:=Order(I);
    //vprintf CRT,2 : "S :=%o;\nI := %o;\nJ := %o;\n\/\/[a,b] =\nelts := %o;\n",PrintSeqAlgEtElt(ZBasis(S)),PrintSeqAlgEtElt(ZBasis(I)),PrintSeqAlgEtElt(ZBasis(J)),PrintSeqAlgEtElt([a,b]); 
    require S eq Order(J): "the ideals must be of the same order";
    require IsCoprime(I,J) : "the ideals must be coprime";
//    I_min:=MinimalInteger(I);
//    J_min:=MinimalInteger(J);
//    g,c1,d1:=XGCD(I_min,J_min);
//    if g ne 1 then
        QIJ,qIJ:=quo_frac_idl(S,I*J);
        QJ,qJ:=quo_frac_idl(S,J);
        QI,qI:=quo_frac_idl(S,I);
        D,pI,pJ,preI,preJ:=DirectSum(QI,QJ);
        zbS:=ZBasis(S);
        X:=[ qIJ(z) : z in zbS];
        Y:=[ pI(qI(z))+pJ(qJ(z)) : z in zbS];
        iso:=Isomorphism(QIJ,D,X,Y);
        qD:=hom<S->D | x:->pI(qI(x))+pJ(qJ(x))>;
//    else
//        //g:=c1*I_min+d1*J_min
//        c:=c1*I_min;
//        d:=d1*J_min;
//        e:=a*d+b*c;
//    end if;
    return qIJ,iso,qD;
end intrinsic;

////////
// CRT
///////

    AttachSpec("~/packages_github/AlgEt/spec");
    import "Ord.m" : crQZ , crZQ , Columns , hnf , MatrixAtoQ , MatrixAtoZ , MatrixQtoA , meet_zbasis , inclusion_matrix;
    _<x>:=PolynomialRing(Integers());
    f:=x^4 + 6*x^2 + 25;
    K:=EtaleAlgebra(f);
    O:=MaximalOrder(K);
    pp:=PrimesAbove(2*O);
    P:=pp[1]; 
    Q:=pp[2];
    mP:=MatrixAtoQ(ZBasis(P));
    mQ:=MatrixAtoQ(ZBasis(Q));
    C:=VerticalJoin(mP,mQ);
    d:=Denominator(C);
    mO:=d*VerticalJoin(MatrixAtoQ(ZBasis(O)),ZeroMatrix(Rationals(),Dimension(K)));
    HH,V:=HermiteForm(crQZ(mO));
    H,U:=HermiteForm(crQZ(d*C)); // U*d*C = H = V*mO
    dK:=d*One(K);

    cc:=Matrix([AbsoluteCoordinates([One(K)],O)[1] cat [0 : i in [1..Dimension(K)]]])*crZQ(V^-1*U); cc;
    assert &+[Eltseq(cc)[i]*(ZBasis(P) cat ZBasis(Q))[i] : i in [1..#Eltseq(cc)]] eq One(K);


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

    /* based on new idea
    // test 3
    t0:=Cputime();
    out3:=[];
    qIJ,iso,qD:=ChineseRemainderTheorem2(pp13[1],pp13[2]);
    for pair in pairs do
        a:=pair[1];
        b:=pair[2];
        e:=(qD(a)+qD(b))@@iso@@qIJ;
        Append(~out3,e);
    end for;
    Cputime(t0);
    assert forall{i : i in [1..#out1] | out1[i] eq out3[i]};
    */


*/
