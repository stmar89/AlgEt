AttachSpec("../../spec");
SetAssertions(2);
time_start:=Cputime();


    printf "### Testing Creation of Algebra:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    SetAssertions(2);
    _<x>:=PolynomialRing(Integers());
    printf ".";
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    assert #Basis(A) eq Dimension(A);

    seq:=[x^2-5,x^2-7];
    seq:=[NumberField(f) : f in seq];
    A:=EtaleAlgebra(seq);
    printf ".";

    seq:=[x-1,x-20];
    A:=EtaleAlgebra(&*seq);
    printf ".";

    SetAssertions(1);
    printf " all good!\n";






    printf "### Testing Attributes and Equality:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    SetAssertions(2);
    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    B:=EtaleAlgebra(f);
    assert not A eq B;
    assert A ne B;
    assert A eq A;
    assert A cmpeq A;
    assert not A cmpeq B;
    SetAssertions(1);
    printf " all good!\n";






    printf "### Testing IsNumberField:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    SetAssertions(2);
    _<x>:=PolynomialRing(Integers());
    printf ".";
    f:=x^8+16;
    A:=EtaleAlgebra(f);
    test,K,i:=IsNumberField(A);
    assert test;
    for j in [1..100] do
        a:=Random(A);
        b:=Random(A);
        assert i(a+b) eq i(a)+i(b);
        assert i(a*b) eq i(a)*i(b);
        ii:=Inverse(i);
        OK:=MaximalOrder(K);
        a:=Random(OK,3);
        b:=Random(OK,3);
        assert ii(a+b) eq ii(a)+ii(b);
        assert ii(a*b) eq ii(a)*ii(b);
    end for;



    seq:=[x^2-5,x^2-7];
    seq:=[NumberField(f) : f in seq];
    A:=EtaleAlgebra(seq);
    assert not IsNumberField(A);
    printf ".";

    SetAssertions(1);
    printf " all good!\n";






    printf "### Testing Elements:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    SetVerbose("AlgEtQElt",2);

    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    a:=PrimitiveElement(A);
    assert MinimalPolynomial(a) eq f;
    _:=A!1;
    _:=A!(1/2);
    _:=Random(A,3)+Random(A,3);
    _:=Random(A,3)-Random(A,3);
    _:=Random(A,3)*Random(A,3);
    _:=Random(A,3)/RandomUnit(A,3);
    _:=Random(A);
    assert not IsIntegral(A!1/2);
    assert IsIntegral(A!2);

    ort:=OrthogonalIdempotents(A);
    idem:=Idempotents(A);
    assert ort subset idem;
    assert One(A) in idem;
    assert Zero(A) in idem;
    assert forall{i : i in idem | i^2 eq i};
    assert forall{i : i,j in [1..#ort] | ort[i]*ort[j] eq ort[i]*KroneckerDelta(i,j)};

    for n in [1..100] do
        a:=Random(A,3);
        coord1:=AbsoluteCoordinates([a],Basis(A));
        coord2:=AbsoluteCoordinates([a],PowerBasis(A));
        printf ".";
    end for;

    seq:=[x^2-5,x^2-7];
    seq:=[NumberField(f) : f in seq];
    A:=EtaleAlgebra(seq);
    a:=PrimitiveElement(A);
    assert MinimalPolynomial(a) eq DefiningPolynomial(A);
    _:=A!1;
    _:=A!(1/2);
    _:=A!<seq[1]!1,seq[2]!(1/2)>;
    _:=A![1,2];
    assert A!1 + A!(1/2) eq A!3/2;
    e:=A!<seq[1]!1,seq[2]!(1/2)>/A![1,2];
    assert not IsIntegral(e);
    for n in [1..100] do
        a:=Random(A,3);
        coord1:=AbsoluteCoordinates([a],Basis(A));
        coord2:=AbsoluteCoordinates([a],PowerBasis(A));
        printf ".";
    end for;

    // testing sequences
    //AttachSpec("~/packages_github/AlgEt/spec");
    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    seq:=[Random(A) : i in [1..10^3]];
    c0:=&+seq;
    c:=seq[1];
    for i in [2..#seq] do
        c+:=seq[i];
    end for;
    assert c eq c0;
    c0:=&*seq;
    c:=seq[1];
    for i in [2..#seq] do
        c*:=seq[i];
    end for;
    assert c eq c0;

    s1:=&+[seq[i]*seq[i] : i in [1..#seq]];
    s2:=SumOfProducts(seq,seq);
    assert s1 eq s2;
    seq0:=[Random(-10,10) : i in [1..#seq]];
    s1:=&+[seq0[i]*seq[i] : i in [1..#seq]];
    s2:=SumOfProducts(seq0,seq);
    s3:=&+[seq[i]*seq0[i] : i in [1..#seq]];
    s4:=SumOfProducts(seq,seq0);
    assert s1 eq s2 and s1 eq s3 and s1 eq s4;
    seq1:=[Random(-10,10)/Random(1,20) : i in [1..#seq]];
    s1:=&+[seq1[i]*seq[i] : i in [1..#seq]];
    s2:=SumOfProducts(seq1,seq);
    s3:=&+[seq[i]*seq1[i] : i in [1..#seq]];
    s4:=SumOfProducts(seq,seq1);
    assert s1 eq s2 and s1 eq s3 and s1 eq s4;

    printf " all good!\n"; 






    printf "### Testing Homs:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    homs:=HomsToC(A);
    a:=PrimitiveElement(A);
    assert &and[ Abs(Evaluate(f,h(a))) lt 10^-20 : h in homs ];
    printf ".";
    old_prec:=Precision(Codomain(homs[1]));
    new_prec:=10*old_prec;
    homs:=HomsToC(A : Prec:=new_prec);
    assert Precision(Codomain(homs[1])) eq new_prec;
    printf ".";
    B:=EtaleAlgebra(Components(A) cat [NumberField(x^2+2)]);
    img:=[ B!(Components(b) cat <0>) : b in AbsoluteBasis(A) ];
    incl:=Hom(A,B,img : CheckMultiplicative:=true );
    assert incl(One(A)) ne One(B);
    printf ".";
    assert forall{ a : a in AbsoluteBasis(A) | MinimalPolynomial(incl(a)) eq MinimalPolynomial(a)};
    printf ".";
    aut:=[ Automorphisms(K) : K in Components(B) ];
    aut:=[ Random(a) : a in aut ];
    img:=[ B!<aut[i](Components(a)[i]) : i in [1..#aut] > : a in AbsoluteBasis(B) ];
    aut:=Hom(B,B,img: CheckMultiplicative:=true, CheckUnital:=true);
    inv:=Inverse(aut);
    assert forall{ b : b in AbsoluteBasis(B) | (inv*aut)(b) eq b and (aut*inv)(b) eq b};
    printf ".";
    A:=EtaleAlgebra(x^2+x+2);
    B:=EtaleAlgebra(x^2-x+2);
    pi:=PrimitiveElement(B);
    m:=Hom(A,B,[(-pi)^i : i in [0..Dimension(B)-1]]);
    assert Inverse(m)(pi) eq -PrimitiveElement(A);
    assert m(One(A)) eq One(B) and Inverse(m)(One(B)) eq One(A);
    printf ".";
    printf " all good!\n"; 






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





    
    printf "### Testing Orders:";
	//AttachSpec("~/packages_github/AlgEt/spec");
    SetVerbose("AlgEtQOrd",1);
    SetAssertions(2);

    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    O1:=Order(Basis(A));
    _:=ZBasis(O1);
    _:=Generators(O1);
    O2:=Order(AbsoluteBasis(A) : Check:=0);
    O3:=Order(AbsoluteBasis(A));
    _:=O1 eq O2;
    _:=O2 eq O3;
    assert EquationOrder(A) ne ProductOfEquationOrders(A);
    printf ".";

    // we test KnownOrders
    OA:=MaximalOrder(A);
    O:=Order(ZBasis(OA));
    assert assigned O`IsMaximal; //at creation O<-OA since we call IsKnownOrder
    assert O eq OA;
    printf ".";

    ff:=Conductor(EquationOrder(A));
    T:=MultiplicatorRing(ff);
    assert T in A`KnownOrders[Index(T)];
    assert assigned T`IsMaximal and assigned T`IsProductOfOrders; //these two attributes are assigned when we created OA above.
                                                                  // hence T points to the same memory spot of OA
    printf ".";


    O:=Order(ZBasis(OA));
    assert IsProductOfOrders(O);
    assert IsMaximal(O);
    printf ".";

    O:=MaximalOrder(A);
    G:=[[ Random(O) : i in [1..3] ] : i in [1..100]];
    S:=[ Order(s) : s in G ];
    _:=#Seqset(S);
    printf ".";

    assert forall{z : z in ZBasis(O1) | z in O1 };
    for O in [O1,O2,O3] do
        for i in [1..100] do
            assert Random(O) in O;
        end for;
    end for;
    printf ".";

    seq:=[x^2-5*25,x^2-7*49];
    seq:=[NumberField(f) : f in seq];
    A:=EtaleAlgebra(seq);
    O1:=Order(Basis(A));
    O2:=Order(AbsoluteBasis(A));
    assert O1 eq O2;
    for O in [O1,O2] do
        for i in [1..100] do
            assert Random(O) in O;
        end for;
    end for;
    printf ".";

    E1:=EquationOrder(A);
    E2:=ProductOfEquationOrders(A);
    E3:=Order(A,<EquationOrder(seq[1]) , MaximalOrder(seq[2])>);
    E4:=Order(A,<MaximalOrder(seq[1]) , EquationOrder(seq[2])>);
    OA:=MaximalOrder(A);
    for E in [E1,E2,E3,E4] do 
        _:=Index(OA,E);
    end for;
    assert E1 subset E2;
    assert E1 subset E3;
    assert E3*E4 eq OA;
    assert OA meet E1 eq E1;
    assert not E2 meet E1 eq E2;
    assert E3 meet E4 eq E2;
    O:=MaximalOrder(A);

    printf " all good!\n"; 






    printf "### Testing Ideals:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    SetVerbose("AlgEtQIdl",1);
    SetAssertions(2);

    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    E1:=EquationOrder(A);
    E2:=ProductOfEquationOrders(A);
    for i in [1..100] do
        a:=Random(A);
        assert 1 eq #{Ideal(E1,a),a*E1,E1*a};
        assert 1 eq #{Ideal(E2,a),a*E2,E2*a};
        assert a*E2 eq E2!!(a*E1);
        assert E1!!(E2!!(a*E1)) ne a*E1;
        printf ".";
    end for;

    _:={ Random(E1)*E1 : i in [1..100] };
    _:=a*E1 + E1!!(a*E2);
    _:=&*[Random(E1)*E1 : i in [1..100]];
    _:=&*[(Random(E1)*E1+Random(E1)*E1) : i in [1..100]];
    l:=[(Random(E1)*E1+Random(E1)*E1) : i in [1..100]];
    assert forall{ I : I in l | not assigned I`ZBasis };
    assert forall{ I : I in l | I eq Ideal(E1,ZBasis(I))};
    assert forall{ I : I in l | assigned I`ZBasis };
    _:=&+[ i eq j select 1 else 0 : i,j in l ];
    I:=Ideal(E1,[Random(E1) : i in [1..100]]);
    J:=&*[I : i in [1..100]];
    JJ:=I^100;
    assert J eq JJ;
    _:=&meet[(Random(E1)*E1+Random(E1)*E1) : i in [1..100]];
    assert forall{ I : I in [Random(E1)*E1 : i in [1..100]] | not IsProductOfIdeals(I)};
    printf ".";
    assert forall{ I : I in [Random(E2)*E2 : i in [1..100]] | IsProductOfIdeals(I)};
    _:=[TraceDualIdeal(Random(E1)*E1+Random(E1)*E1) : i in [1..100]];
    _:=[IsIntegral(Random(E1)*E1+Random(E1)*E1) : i in [1..100]];
    _:=[MakeIntegral(Random(E1)*E1+Random(E1)*E1) : i in [1..100]];
    printf ".";

    ids:=[ Ideal(E1,[Random(E1) : i in [1..10]]) : i in [1..200]]; 
    _:=&meet(ids);
    rr:=[ResidueRing(E1,I) : I in ids ];
    cc:=[ ColonIdeal(I,J) : I,J in ids[1..10]  ];

    O:=MaximalOrder(A);
    test,O_prod:=IsProductOfOrders(O);
    assert test;
    for i in [1..100] do
        _:=Ideal(O,<ideal< O_prod[i] | [Random(O_prod[i],3): j in [1..3]]> : i in [1..#O_prod] >);
        printf ".";
    end for;

    ids:=[ Ideal(E1,[Random(E1) : i in [1..10]]) : i in [1..20]]; 
    cc2:=[ TraceDualIdeal(TraceDualIdeal(I)*J) : I,J in ids  ];
    cc:=[ ColonIdeal(I,J) : I,J in ids  ];
    assert cc eq cc2;
    printf ".";

    // testing KnownPowers 
    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    E1:=EquationOrder(A);
    I:=Conductor(E1);
    assert not IsInvertible(I);
    _:=I^2;
    is_def,I2:=IsDefined(I`KnownPowers,2);
    assert is_def and I2 eq I*I;
    _:=I^5;
    is_def,I4:=IsDefined(I`KnownPowers,4);
    assert is_def and I4 eq (I*I)*(I*I);
    is_def,I5:=IsDefined(I`KnownPowers,5);
    assert is_def and I5 eq (I*I)*(I*I)*I;

    SetAssertions(1);
    printf " all good!\n";







  printf "### Testing ZBasisLLL:";
	//AttachSpec("~/packages_github/AlgEt/spec");
	SetAssertions(2);
	_<x>:=PolynomialRing(Integers());
  f:=x^4-100*x^3-100*x^2-100*x-100;
  K:=EtaleAlgebra(f);
  E:=EquationOrder(K);
  pp:=SingularPrimes(E);
  I:=&*(pp);
  J:=&*(pp);
  ZBasisLLL(I);
  assert ZBasis(J) ne ZBasis(I);
  assert J eq I;
  printf " all good!\n"; 






    printf "### Testing Quotients:";
    //AttachSpec("~/packages_github/AlgEt/spec");
	SetAssertions(2);
	_<x>:=PolynomialRing(Integers());
    f:=(x^4+16);
	A:=EtaleAlgebra(f);
	E:=EquationOrder(A);
    icm:=ICM(E);
    pp:=SingularPrimes(E);
    for P in pp do
        r:=ResidueField(P);
        _:=PrimitiveElementResidueField(P);
        assert #r eq #Quotient(OneIdeal(E),P);
        assert #r eq #ResidueRing(E,P);
        assert #r eq #QuotientVS(OneIdeal(E),P,P);
        printf ".";
    end for;

    for I in icm do
        _:=Quotient(OneIdeal(E),MakeIntegral(I));
        for P in pp do
            PI:=P*I;
            n:=#Quotient(I,PI);
            assert n eq #QuotientVS(I,PI,P);
            printf ".";
        end for;
    end for;

    printf " all good!\n"; 






    printf "### Testing OverOrders:";
	//AttachSpec("~/packages_github/AlgEt/spec");

    SetVerbose("OverOrders",1);
    SetAssertions(1);

    _<x>:=PolynomialRing(Integers());
    f:=(x^4+16);
    A:=EtaleAlgebra(f);
    O:=MaximalOrder(A);
    assert #FindOverOrders(O) eq 1;
    assert #MinimalOverOrders(O) eq 0;

    E:=EquationOrder(A);
    oo:=FindOverOrders(E);
    assert #oo eq 11;

    printf " all good!\n"; 





    
    printf "### Testing GraphOverOrders:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    P<x>:=PolynomialRing(Integers());
    fs:=[ 
          x^8 + 16, //1 sing prime
          x^4-10000*x^3-10000*x^2-10000*x-10000, //2 sing primes
          x^4-30^3*x^3-30^3*x^2-30^3*x-30^3 //3 sing primes
        ];
    for f in fs do
        printf ".";
        A:=EtaleAlgebra(f);
        _:=GraphOverOrders(MaximalOrder(A));
        R:=EquationOrder(A); 
        _:=OverOrders(R);
        SetAssertions(2);
        _:=GraphOverOrders(R);
        SetAssertions(1);
    end for;
    printf " all good!\n";






    printf "### Testing Trace and Norm:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    SetVerbose("AlgEtQTraceNorm",1);

    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    for i in [1..100] do
        a:=Random(A);
        b:=Random(A);
        assert Trace(a)+Trace(b) eq Trace(a+b);
        assert Norm(a)*Norm(b) eq Norm(a*b);
    end for;
    printf " all good!\n"; 





    
    printf "### Testing Completion:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    PP<x>:=PolynomialRing(Integers());
    polys:=[
        x^6+3*x^4-10*x^3+15*x^2+125,
        (x^2+5)*(x^4-4*x^3+5*x^2-20*x+25),
        (x^4-5*x^3+15*x^2-25*x+25)*(x^4+5*x^3+15*x^2+25*x+25)
        ];
    for h in polys do
        L:=EtaleAlgebra(h);
        a:=PrimitiveElement(L);
        O:=MaximalOrder(L);
        p:=5;
        pp:=PrimesAbove(p*O);
        for P in pp do
            C,mC:=Completion(P);
        end for;
        printf ".";
    end for;
    printf " all good!\n";






    printf "### Testing Complex Conjugation:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    _<x>:=PolynomialRing(Integers());
    f:=x^4 + 6*x^2 + 25;
    K:=EtaleAlgebra(f);
    pi:=PrimitiveElement(K);
    pib:=ComplexConjugate(pi);
    assert pib eq 5/pi;
    R:=Order([pi,pib]);
    oo:=FindOverOrders(R);
    O:=MaximalOrder(K);
    S:=[ S : S in oo | Index(O,S) eq 8 ][1]; //there is only one order of index 8, so conjugate stable
    assert IsConjugateStable(R);
    assert IsConjugateStable(TraceDualIdeal(R));
    assert IsConjugateStable(O);
    assert IsConjugateStable(TraceDualIdeal(O));
    assert IsConjugateStable(S);
    assert IsConjugateStable(TraceDualIdeal(S));
    assert not IsConjugateStable(EquationOrder(K));
    printf ".";
    printf " all good!\n";






    printf "### Testing CM-types:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    _<x>:=PolynomialRing(Integers());
    polys:=[
    x^4+x^2+529,
    x^4+11*x^3+73*x^2+319*x+841,
    x^4-4*x^3+8*x^2-116*x+841,
    x^4+4*x^3+8*x^2+116*x+841,
    x^4-17*x^2+841,
    x^4-x^3+26*x^2-31*x+961,
    x^4-6*x^3+43*x^2-186*x+961,
    x^4-4*x^3+8*x^2-124*x+961,
    x^4+2*x^3+52*x^2+62*x+961,
    x^4+x^3+26*x^2+31*x+961
    ];

    for f in polys do
        A:=EtaleAlgebra(f);
        all:=AllCMTypes(A);
        _:=[ CMPositiveElement(PHI) : PHI in all ];
        for i,j in [1..#all] do
            assert (i eq j) eq (all[i] eq all[j]);
        end for;
        for i,j in [1..#all] do
            assert (i eq j) eq (all[i] eq ChangePrecision(all[j],2*Precision(all[j])));
        end for; 
        ChangePrecision(~all[1],60);
        assert Precision(all[1]) eq 60;
        printf ".";
    end for;
    printf " all good!\n";






    printf "### Testing IntermediateIdeals:";
	SetAssertions(2);
	_<x>:=PolynomialRing(Integers());
    f:=x^3-100*x^2-100*x-100;
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    O:=MaximalOrder(K);
    ff:=Conductor(E);
    _:=MinimalIntermediateIdeals(E!!OneIdeal(O),ff);
    printf ".";
    _:=IntermediateIdeals(E!!OneIdeal(O),ff);
    printf ".";
    _:=IntermediateIdealsWithPrescribedMultiplicatorRing(E!!OneIdeal(O),ff);
    printf ".";
    _:=MaximalIntermediateIdeals(E!!OneIdeal(O),ff);
    printf ".";
    _:=IntermediateIdealsWithTrivialExtension(E!!OneIdeal(O),ff,O);
    printf ".";
    _:=IntermediateIdealsWithTrivialExtensionAndPrescribedMultiplicatorRing(E!!OneIdeal(O),ff,O);
    printf ".";
    f:=x^4+291*x^3-988*x^2-1000*x-1000;
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    O:=MaximalOrder(K);
    ff:=Conductor(E);
    _:=MinimalIntermediateIdeals(E!!OneIdeal(O),ff);
    printf ".";
    _:=MaximalIntermediateIdeals(E!!OneIdeal(O),ff);    
    printf ".";
    _:=IntermediateIdealsWithTrivialExtension(E!!OneIdeal(O),ff,O);
    printf ".";
    _:=IntermediateIdealsWithTrivialExtensionAndPrescribedMultiplicatorRing(E!!OneIdeal(O),ff,O);
    printf ".";
    SetAssertions(1);
    printf " all good!\n"; 






    printf "### Testing IdealsOfIndex:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    SetVerbose("IdealsOfIndex",1);
    SetAssertions(2);

    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    E:=EquationOrder(A);
    O:=MaximalOrder(A);
    ind:=Index(O,E);
    for N in [1..15] do
        printf "."; 
        // test with maximal order
        assert Seqset(IdealsOfIndex(O,N)) eq Seqset(IdealsOfIndex(O,N : Method:="Slow"));

        //test with equation order
        if IsCoprime(N,ind) then
            assert Seqset(IdealsOfIndex(E,N)) eq Seqset(IdealsOfIndex(E,N : Method:="Slow"));
        else
            _:=IdealsOfIndex(E,N);
        end if;
    end for;
    SetAssertions(1);    
    printf " all good!\n"; 





    
    printf "### Testing ShortEltSmallRep:";
	//AttachSpec("~/packages_github/AlgEt/spec");
	_<x>:=PolynomialRing(Integers());
    f:=(x^2+5)*(x^2+7)*(x^2+11);
    assert IsSquarefree(f);
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    ff:=Conductor(E);
    _:=ShortElement(ff);
    oo:=FindOverOrders(E); 
    for S in oo do
        printf ".";
        ff:=Conductor(S);
        _:=ShortElement(ff);
    end for;
    printf " all good!\n"; 






    printf "### Testing MinimalGenerators:";
	//AttachSpec("~/packages_github/AlgEt/spec");
    SetClassGroupBounds("GRH");
	_<x>:=PolynomialRing(Integers());
    f:=x^4-1000*x^3-1000*x^2-1000*x-1000;
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    P,p:=PicardGroup(E : GRH:=true); //~10 secs

    for g in Generators(P) do 
        I:=p(g);
        TwoGeneratingSet(I);
        assert #Generators(I) le 2;
        printf ".";
    end for;


    // test if TwoGeneratingSet makes the power faster
    // Conlcusion: yes. By quite a bit!
    f:=x^4-100*x^3-100*x^2-100*x-100;
    A:=EtaleAlgebra(f);
	E:=EquationOrder(A);
    P,p:=PicardGroup(E : GRH:=true);
    repeat
        Ii:=Random(P);
    until Ii ne Zero(P);
    I:=p(Ii);

    delete I`IsInvertible;
    exp:=[ Random(2,30) : i in [1..100]];
    l1:=[ I^i : i in exp ];
    printf ".";

    assert IsInvertible(I);
    TwoGeneratingSet(I);
    assert #Generators(I) eq 2;
    l2:=[ I^i : i in exp ];
    assert l1 eq l2;
    printf ".";

    I:=SmallRepresentative(I);
    delete I`IsInvertible;
    l1:=[ I^i : i in exp ];

    assert IsInvertible(I);
    TwoGeneratingSet(I);
    l2:=[ I^i : i in exp ];
    assert l1 eq l2;
    printf " all good!\n"; 






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






    printf "### Testing PicardGroup and UnitGroup:";
    //AttachSpec("~/packages_github/AlgEt/spec");
	SetAssertions(2);
	_<x>:=PolynomialRing(Integers());

	// very fast
    f:=(x^4+16);
	A:=EtaleAlgebra(f);
	E:=EquationOrder(A);
	oo:=FindOverOrders(E);
    for S in oo do
        P,p:=PicardGroup(S);
        U,u:=UnitGroup(S);
        printf ".";	
    end for;

	//AttachSpec("~/packages_github/AlgEt/spec");
    _<x>:=PolynomialRing(Integers());
    f:=x^4-1000*x^3-1000*x^2-1000*x-1000;
    SetAssertions(1);
    SetVerbose("AlgEtQPicardGroup",1);
    SetVerbose("AlgEtQIdl",1);
    SetVerbose("ShortEltSmallRep",1);
    A:=EtaleAlgebra(f);
    E:=EquationOrder(A);
    P,p:=PicardGroup(E : GRH:=true);
    U,u:=UnitGroup(E : GRH:=true);
    assert #P eq 3548000;
    printf ".";
    SetAssertions(1);

	//AttachSpec("~/packages_github/AlgEt/spec");
    _<x>:=PolynomialRing(Integers());
    f:=x^8+16;
    A:=EtaleAlgebra(f);
    OA:=MaximalOrder(A);
    P:=PrimesAbove(2*OA)[1];
    E:=Order(ZBasis(P^2));
    assert not IsMaximal(E);
    P:=Conductor(E);
    assert IsPrime(P);
    assert #ResidueRingUnits(P) eq Index(E,P)-1;
    assert forall{ i : i in [1..20] | #ResidueRingUnits(Pi) eq Index(E,Pi) - Index(P,Pi) where Pi:=P^i};

    printf " all good!\n"; 






    printf "### Testing Primes and Factorizaton:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    SetAssertions(2);

    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16);
    A:=EtaleAlgebra(f);
    E1:=EquationOrder(A);
    ff:=Conductor(E1);
    assert PrimesAbove(Conductor(E1)) eq SingularPrimes(E1);
    printf ".";

    f:=x^6 + 8*x^5 + 50*x^4 + 200*x^3 + 1250*x^2 + 5000*x + 15625;
    A:=EtaleAlgebra(f);
    assert #SingularPrimes(MaximalOrder(A)) eq 0;
    R:=EquationOrder(A);
    assert #SingularPrimes(R) eq 5;
    assert #NonInvertiblePrimes(R) eq 5;
    printf ".";

    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    E1:=EquationOrder(A);
    E2:=ProductOfEquationOrders(A);
    
    _:=PrimesAbove(Conductor(E1));
    _:=PrimesAbove(Conductor(E2));
    assert IsGorenstein(E1);
    assert IsGorenstein(E2);

    ids:=[ Ideal(E1,[Random(E1) : i in [1..10]]) : i in [1..100]];
    ids0:=[ I : I in ids | I ne OneIdeal(E1) and IsInvertible(I) ];
    ids:=[];
    for I in ids0 do
        _,J:=CoprimeRepresentative(I,Conductor(E1));
        Append(~ids,J);
    end for;
    facs:=[ Factorization(I) : I in ids ];
    printf ".";
    SetAssertions(1);
    printf " all good!\n"; 










    printf "### Testing LowCohenMacaulayType:";
    SetAssertions(2);
    _<x>:=PolynomialRing(Integers());
    f:=x^4-10000*x^3-10000*x^2-10000*x-10000; 
    //AttachSpec("~/packages_github/AlgEt/spec");
    A:=EtaleAlgebra(f);
    E:=EquationOrder(A);
    oo:=FindOverOrders(E); // ~13 secs
    assert #oo eq 297;
    for S in oo do
        pp:=SingularPrimes(S);
        for P in pp do        
            _:=CohenMacaulayTypeAtPrime(S,P);
        end for;
        pp_ng:=NonGorensteinPrimes(S);
        assert pp_ng subset pp;
        _:=CohenMacaulayType(S);
        printf ".";
    end for;
    SetAssertions(1);
    printf " all good!\n"; 






    printf "### Testing WKICM:";
	//AttachSpec("~/packages_github/AlgEt/spec");
	SetAssertions(2);
	_<x>:=PolynomialRing(Integers());
    f:=x^4-1000*x^3-1000*x^2-1000*x-1000;
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    assert #WKICM(E) eq 25;

    f:=x^4+291*x^3-988*x^2-1000*x-1000;
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    assert #WKICM(E) eq 20;

    f:=x^3+31*x^2+43*x+77;
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    assert #FindOverOrders(E) eq 15;
    assert #WKICM(E) eq 23;
    SetAssertions(1);
    printf " all good!\n"; 






    printf "### Testing WKEq:";
	//AttachSpec("~/packages_github/AlgEt/spec");
	SetAssertions(2);
	_<x>:=PolynomialRing(Integers());
    f:=x^4-100*x^3-100*x^2-100*x-100;
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    assert not IsWeakEquivalent(E,Conductor(E));
    assert not IsWeakEquivalent(OneIdeal(E),Conductor(E));
    assert not IsWeakEquivalent(E,MaximalOrder(K));
    assert IsWeakEquivalent(OneIdeal(MaximalOrder(K)),Conductor(E));
    SetAssertions(1);
    printf " all good!\n"; 
    






    printf "### Testing ICM:";
	SetAssertions(2);
    SetClassGroupBounds("GRH");
	_<x>:=PolynomialRing(Integers());
    polys:=[
        x^3-100*x^2-100*x-100,
        x^4+291*x^3-988*x^2-1000*x-1000,
        x^3+31*x^2+43*x+77
        ];
    for f in polys do
        K:=EtaleAlgebra(f);
        E:=EquationOrder(K);
        _:=ICM(E);
    end for;
    SetAssertions(1);
    printf " all good!\n";






    printf "### Testing TotRealPos:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    PP<x>:=PolynomialRing(Integers());
    SetAssertions(2);

    f:=x^8+16;
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    R:=Order([F,2/F]);
    oo:=FindOverOrders(R);
    for iS in [1..#oo] do
        S:=oo[iS];
        _:=TotallyRealUnitGroup(S);
        _:=TotallyRealPositiveUnitGroup(S);
    end for;
    SetAssertions(1);
    printf " all good!\n"; 





    
    printf "### Testing Print Saving:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    _<x>:=PolynomialRing(Integers());
    f:=x^6 - 3*x^5 - 3*x^4 + 65*x^3 - 48*x^2 - 768*x + 4096;
    A:=EtaleAlgebra(f);
    E:=EquationOrder(A);
    seq,str:=PrintSeqAlgEtQElt(ZBasis(E));
    assert Order([ A! s : s in eval(str)]) eq E;
    printf ".";

    //AttachSpec("~/packages_github/AlgEt/spec");
    _<x>:=PolynomialRing(Integers());
    f:=x^8+16;
    A:=EtaleAlgebra(f);
    O:=MaximalOrder(A);
    str:=PrintWKICM(O);
    O1:=LoadWKICM(str);
    printf ".";

    //AttachSpec("~/packages_github/AlgEt/spec");
    _<x>:=PolynomialRing(Integers());
    f:=x^8+16;
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    R:=Order([F,2/F]);
    str:=PrintWKICM(R);
    R1:=LoadWKICM(str);
    assert #WKICM(R) eq #WKICM(R1);
    assert #FindOverOrders(R) eq #FindOverOrders(R1);
    printf ".";

    printf " all good!\n"; 




Cputime(time_start);
quit;