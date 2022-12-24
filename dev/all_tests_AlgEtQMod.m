/* vim: set syntax=magma :*/
    "-------------------------------------------------------------";
    "-------------------------------------------------------------";
    "All tests for AlgEtQMod";
    time_start_slow_tests:=Cputime();
    
    "-------------------------------------------------------------";
    "-------------------------------------------------------------";
    "--------------------------Modules.m--------------------------";
    "-------------------------------------------------------------";
    "-------------------------------------------------------------";

    AttachSpec("~/packages_github/AlgEt/spec");
    AttachSpec("~/packages_github/AlgEt/specMod");
    _<x>:=PolynomialRing(Integers());
    m1:=x^4 - 2*x^2 + 9;
    m2:=x^2 -5*x + 7;
    K1:=NumberField(m1);
    K2:=NumberField(m2);
    K:=EtaleAlgebra([K1,K2]);
    V:=EtaleAlgebra([K1,K2,K1,K2,K2]);
    m:=NaturalAction(K,V);
    F:=PrimitiveElement(K);
    [MinimalPolynomial(c) : c in Components(m(F))];

    Vnf,Ve:=Components(V);
    gens:=&cat[ [Ve[i](z) : z in Basis(MaximalOrder(Vnf[i])) ] : i in [1..#Vnf]];
    gens;
    O:=MaximalOrder(K);
    M:=Module(O,m,gens);
    N:=Module(O,m,< 1*MaximalOrder(Vnf[i]) : i in [1..#Vnf] >);
    Generators(M);
    ZBasis(M) eq ZBasis(N);
    assert N eq M;
    E:=EquationOrder(K);
    assert not IsMaximal(E);
    NE:=Module(E,m,< 1*MaximalOrder(Vnf[i]) : i in [1..#Vnf] >);
    assert N subset NE and NE subset N; // N has multiplicatorring O
    N2:=Module(E,m,< 4*MaximalOrder(Vnf[i]) : i in [1..#Vnf] >);
    assert N2 ne NE;
    assert N2 subset NE;
    assert not NE subset N2;
    _:=Quotient(N,NE); 

    "-------------------------------------------------------------";
    "-------------------------------------------------------------";
    "-------------------IntermediateModules.m----------------------";
    "-------------------------------------------------------------";
    "-------------------------------------------------------------";

    AttachSpec("~/packages_github/AlgEt/spec");
    AttachSpec("~/packages_github/AlgEt/specMod");
    _<x>:=PolynomialRing(Integers());
    m1:=x^4 - 2*x^2 + 9;
    m2:=x^2 -5*x + 7;
    K1:=NumberField(m1);
    K2:=NumberField(m2);
    K:=EtaleAlgebra([K1,K2]);
    V:=EtaleAlgebra([K1,K2,K2]);
    m:=NaturalAction(K,V);
    F:=PrimitiveElement(K);
    [MinimalPolynomial(c) : c in Components(m(F))];

    Vnf,Ve:=Components(V);
    gens:=&cat[ [Ve[i](z) : z in Basis(MaximalOrder(Vnf[i])) ] : i in [1..#Vnf]];
    gens;
    O:=MaximalOrder(K);
    M:=Module(O,m,gens);
    N:=Module(O,m,< 1*MaximalOrder(Vnf[i]) : i in [1..#Vnf] >);
    Generators(M);
    ZBasis(M) eq ZBasis(N);
    assert N eq M;
    E:=EquationOrder(K);
    assert not IsMaximal(E);
    NE:=Module(E,m,< 1*MaximalOrder(Vnf[i]) : i in [1..#Vnf] >);
    assert N subset NE and NE subset N; // N has multiplicatorring O
    N2:=Module(E,m,< 2*MaximalOrder(Vnf[i]) : i in [1..#Vnf] >);
    assert N2 ne NE;
    assert N2 subset NE;
    assert not NE subset N2;
    Q,q:=Quotient(NE,N2);
    #Q;
    time #MinimalIntermediateModules(NE,N2); // 0.1 secs
    time l:=IntermediateModules(NE,N2); #l; // 5 secs
    time l3:=IntermediateModulesWithTrivialExtension(NE,N2,O); // ~1 secs
    NEO:=O!!NE;
    assert #l3 eq #[ M : M in l | O!!M eq NEO ];

    // a much bigger test!
    AttachSpec("~/packages_github/AlgEt/spec");
    AttachSpec("~/packages_github/AlgEt/specMod");
    _<x>:=PolynomialRing(Integers()); 
    h:=x^8 - 4*x^6 + 22*x^4 - 36*x^2 + 81;
    q:=Integers() ! Truncate( ConstantCoefficient(h)^(2/Degree(h)) );
    fac:=Factorization(h);
    g:=&*[f[1] : f in fac];
    K:=EtaleAlgebra(g);
    F:=PrimitiveElement(K);
    V:=q/F;
    ZFV:=Order([F,V]);
    oo:=FindOverOrders(ZFV);
    // The situation is the following: (we denote inclusion by c)
    // ZFV c S c T c O,
    // each inclusion with index 2.
    // S has CohenMacaulayType 2, while ZFV,T and O are Gorenstein.
    // All have only one prime above 2.
    // All have trivial PicardGroup
    // If ZFV were Bass then we woulc have the followin isomorphism classes of AVs in IsogenClass(h):
    ss:={1,2,3,4};
    #[ c : c in car<ss,ss,ss,ss> | c[1] le c[2] and c[2] le c[3] and c[3] le c[4] ];
    // we get 35 classes.
    // But ZFV is NOT Bass!/
    assert #PicardGroup(ZFV) eq 1;

    R:=ZFV;
    O:=MaximalOrder(K);
    nf:=Components(K);
    V:=EtaleAlgebra(&cat[nf : i in [1..2]]);
    Vnf:=Components(V);
    m:=NaturalAction(K,V);
    MO:=Module(R,m,<1*MaximalOrder(Vnf[i]) : i in [1..#Vnf]>);
    MOO:=O!!MO;
    ff:=O!!Conductor(R);
    test,ff_prod:=IsProductOfIdeals(ff);
    assert test;
    // the following works for this exmaple, not in general.
    Mff:=Module(R,m,<ff_prod[1] : i in [1..2]>);
    gensMff:=Generators(Mff);
    time l0:=IntermediateModulesWithTrivialExtension(MO,Mff,O);
    time l1:=IntermediateModules(MO,Mff);
    time l2:=[ I : I in l1 | O!!I eq MO ];
    assert #l2 eq #l0;
    "then do the julia thing...";

    "-------------------------------------------------------------";
    "-------------------------------------------------------------";
    "-----------------------IsomModules.m-------------------------";
    "-------------------------------------------------------------";
    "-------------------------------------------------------------";
    
    // ##########
    // Fast Tests
    // ##########

    AttachSpec("~/packages_github/AlgEt/spec");
    AttachSpec("~/packages_github/AlgEt/specMod");
    SetVerbose("IsomModules",1);
    _<x>:=PolynomialRing(Integers()); 
    g:=x^2+5;
    K:=EtaleAlgebra(g);
    nf:=Components(K);
    pi:=PrimitiveElement(K);
    R:=Order([7*pi]);
    m:=NaturalAction(K,K);
    time icm:=ICM(R);
    time classes:=IsomorphismClasses(R,m : Method:="Magma");
    assert #classes eq #icm; // Since R is Bass, the ICM is Pic(R) cup Pic(O)
    time classes:=IsomorphismClasses(R,m : Method:="julia -J /tmp/Hecke.so ~/packages_github/AlgEt/AlgEtQMod/"); // this is the ICM
    assert #classes eq #icm; // Since R is Bass, the ICM is Pic(R) cup Pic(O)

    // Pic(O) is non trivial
    AttachSpec("~/packages_github/AlgEt/spec");
    AttachSpec("~/packages_github/AlgEt/specMod");
    SetVerbose("IsomModules",1);
    _<x>:=PolynomialRing(Integers()); 
    g:=x^2+15;
    K:=EtaleAlgebra(g);
    nf:=Components(K);
    pi:=PrimitiveElement(K);
    R:=Order([pi]);
    m:=NaturalAction(K,K);
    time classes:=IsomorphismClasses(R,m : Method:="Magma");
    time classes:=IsomorphismClasses(R,m : Method:="julia -J /tmp/Hecke.so ~/packages_github/AlgEt/AlgEtQMod/"); // this is the ICM
    time icm:=ICM(R);
    assert #classes eq #icm; // Since R is Bass, the ICM is Pic(R) cup Pic(O)
    // V = K^2
    V:=EtaleAlgebra(&cat[nf : i in [1..2]]);
    m:=NaturalAction(K,V);
    time classes:=IsomorphismClasses(R,m : Method:="Magma");
    assert #classes eq 6; // because R is Bass 
    time classes:=IsomorphismClasses(R,m : Method:="julia -J /tmp/Hecke.so ~/packages_github/AlgEt/AlgEtQMod/");
    assert #classes eq 6; // because R is Bass 

    AttachSpec("~/packages_github/AlgEt/spec");
    AttachSpec("~/packages_github/AlgEt/specMod");
    SetVerbose("IsomModules",1);
    _<x>:=PolynomialRing(Integers());
    SetVerbose("IsomModules",2);
    m1:=x^2-x+3;
    m2:=x^2+x+3;
    s1:=2;
    s2:=1;
    h:=m1^s1*m2^s2; h;
    K1:=NumberField(m1);
    K2:=NumberField(m2);
    K:=EtaleAlgebra([K1,K2]); // K = K1 x K2
    pi:=PrimitiveElement(K); 
    R:=Order([pi]);
    O:=MaximalOrder(K);
    V:=EtaleAlgebra([K1,K1,K2]); // V = K1^s1 x K2^s2
    m:=NaturalAction(K,V); // m:K -> V component-wise diagonal action of K on V
    time classes:=IsomorphismClasses(R,m : Method:="julia -J /tmp/Hecke.so ~/packages_github/AlgEt/AlgEtQMod/"); // changes this line accordingly to wheter you have used Hecke.Build() or not,
                                                                                                       // and to the appriopriate path to the the packages AlgEtQ
    assert #classes eq 4;

    AttachSpec("~/packages_github/AlgEt/spec");
    AttachSpec("~/packages_github/AlgEt/specMod");
    SetVerbose("IsomModules",1);
    _<x>:=PolynomialRing(Integers()); 
    g:=x^6-x^5+2*x^4-2*x^3+4*x^2-4*x+8;
    K:=EtaleAlgebra(g);
    nf:=Components(K);
    q:=2;
    pi:=PrimitiveElement(K);
    R:=Order([pi,q/pi]);
    O:=MaximalOrder(K);
    assert IsBass(R);
    m:=NaturalAction(K,K);
    time classes:=IsomorphismClasses(R,m : Method:="julia -J /tmp/Hecke.so ~/packages_github/AlgEt/AlgEtQMod/"); // this is the ICM
    time icm:=ICM(R);
    assert #classes eq #icm; // Since R is Bass, the ICM is Pic(R) cup Pic(O)

    // ########## 
    // Slow Tests
    // #########

    // Pic(O) is non triviali. Class construction is rther fast, isomorphism sieveng requires ~10000 secs using julia.
    AttachSpec("~/packages_github/AlgEt/spec");
    AttachSpec("~/packages_github/AlgEt/specMod");
    SetVerbose("IsomModules",1);
    _<x>:=PolynomialRing(Integers()); 
    g:=x^2+15;
    K:=EtaleAlgebra(g);
    nf:=Components(K);
    pi:=PrimitiveElement(K);
    R:=Order([pi]);
    // V = K^3 this is expected to be considerably slower. Because of the bugs described below, I never got to the end.
    // There is a bug also in the Isomorphism testing routing the Magma. It returns false, when it should be true. It has been reported, but I don't expect any quick action.
    V:=EtaleAlgebra(&cat[nf : i in [1..3]]);
    m:=NaturalAction(K,V);
    //time classes:=IsomorphismClasses(R,m : Method:="Magma"); // 800 secs. Faster than expected.
    //assert #classes eq 8; // because R is Bass. Note that the Magma method is bugged.
    time classes:=IsomorphismClasses(R,m : Method:="julia -J /tmp/Hecke.so ~/packages_github/AlgEt/AlgEtQMod/");
    assert #classes eq 8; // because R is Bass 

    AttachSpec("~/packages_github/AlgEt/spec");
    AttachSpec("~/packages_github/AlgEt/specMod");
    SetVerbose("IsomModules",1);
    _<x>:=PolynomialRing(Integers()); 
    g:=x^6-x^5+2*x^4-2*x^3+4*x^2-4*x+8;
    K:=EtaleAlgebra(g);
    nf:=Components(K);
    q:=2;
    pi:=PrimitiveElement(K);
    R:=Order([pi,q/pi]);
    // never finished so far. There are ~1770 candidates.
    // V = K^3
    V:=EtaleAlgebra(&cat[nf : i in [1..3]]);
    m:=NaturalAction(K,V);
    time classes:=IsomorphismClasses(R,m : Method:="julia -J /tmp/Hecke.so ~/packages_github/AlgEt/AlgEtQMod/");
    assert #classes eq 6; // Example 6.1 in "Computing abelian varieties over finite fields isogenous to a power", by Marseglia