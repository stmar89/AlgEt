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

    // fast tests from remote
    screen -r fast_tests;

    quit;
    cd ~/packages_github/AlgEt/dev/fast_tests_AlgEtQ_make/
    git pull; sleep 1;
    make;
    
    // slow tests from remote
    screen -S slow_tests;
    screen -r slow_tests;

    quit;
    cd ~/packages_github/AlgEt/dev/
    git pull; sleep 1;
    magma -b slow_tests_AlgEtQ.m

    // examples AlgEtQ
    screen -r examples_papers_AlgEtQ;

    quit;
    cd ~/packages_github/AlgEt/examples/
    git pull; sleep 1;
    magma -b ideal_class_monoid.txt

    // examples Modules
    screen -S examples_Modules;
    screen -r examples_Modules;

    quit;
    cd ~/packages_github/AlgEt/
    git pull; sleep 1;
    magma -b ~/packages_github/AlgEt/examples/modules_conjugacy_AVs.txt
    magma -b ~/packages_github/AlgEt/dev/all_tests_AlgEtQMod.m



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

    // wkimc _better
    // as a recursion with early exit if type le 2
    wkicm_bar:=function(S)
        ooS:=[ T : T in FindOverOrders(S) | T ne S and T eq MultiplicatorRing(T!!St)];
        ind:=[ Index(T,S) : T in ooS];
        min,pos:=Min(ind);
        assert #[ i : i in ind | i eq min ] eq 1;
        T:=ooS[pos];
        wkT:=[ S!! I : I in WKICM_bar(T) ];
        ff:=ColonIdeal(OneIdeal(S),S!!OneIdeal(T));
        l:=&join{@ IntermediateIdealsWithTrivialExtensionAndPrescribedMultiplicatorRing(I,ff*I,T) : I in wkT @};
        ww:={@ @};
        for I in l do 
            if not exists{J : J in ww | IsWeakEquivalent(I,J)} then
                Include(~ww,I); 
            end if; 
        end for;
        return ww;
    end function;

*/
