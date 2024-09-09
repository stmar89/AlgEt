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

    SetDebugOnError(true);

    path_to_package:="~/"; //folder where the package is installed.

    "-------------------------------------------------------------";
    "-------------------------------------------------------------";
    "All tests for AlgEtQMod";
    time_start_slow_tests:=Cputime();
    
    "-------------------------------------------------------------";
    "-------------------------------------------------------------";
    "--------------------------Modules.m--------------------------";
    "-------------------------------------------------------------";
    "-------------------------------------------------------------";

    AttachSpec(path_to_package cat "AlgEt/spec");
    AttachSpec(path_to_package cat "AlgEt/specMod");
    AttachSpec(path_to_package cat "AlgEt/specMtrx");
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
    _:=Generators(M);
    assert ZBasis(M) eq ZBasis(N);
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
    "-------------------IntermediateModules.m---------------------";
    "-------------------------------------------------------------";
    "-------------------------------------------------------------";

    AttachSpec(path_to_package cat "AlgEt/spec");
    AttachSpec(path_to_package cat "AlgEt/specMod");
    AttachSpec(path_to_package cat "AlgEt/specMtrx");
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

    "-------------------------------------------------------------";
    "-------------------------------------------------------------";
    "------------------IsomModules.m Fast Tests-------------------";
    "-------------------------------------------------------------";
    "-------------------------------------------------------------";
    
    // ##########
    // Fast Tests
    // ##########

    "A very quick test.";
    AttachSpec(path_to_package cat "AlgEt/spec");
    AttachSpec(path_to_package cat "AlgEt/specMod");
    AttachSpec(path_to_package cat "AlgEt/specMtrx");
    SetVerbose("IsomModules",1);
    _<x>:=PolynomialRing(Integers()); 
    g:=x^2+5;
    K:=EtaleAlgebra(g);
    nf:=Components(K);
    pi:=PrimitiveElement(K);
    R:=Order([7*pi]);
    m:=NaturalAction(K,K);
    time icm0:=ICM(R);
    time icm:=IsomorphismClasses(R,m : Method:="Magma", UseSpecializedMethod:=true);
    assert #icm eq #icm0;
    time classes:=IsomorphismClasses(R,m : Method:="Magma", UseSpecializedMethod:=false);
    assert #classes eq #icm; // Since R is Bass, the ICM is Pic(R) cup Pic(O)
    assert forall{i:i,j in [1..#classes]|(i eq j) eq IsIsomorphic(classes[i],classes[j] : UseSpecializedMethod:=true)};
    time classes:=IsomorphismClasses(R,m : Method:="julia -J /tmp/Hecke.so " cat path_to_package cat "AlgEt/AlgEtQMod/", UseSpecializedMethod:=false); // this is the ICM
    assert #classes eq #icm; // Since R is Bass, the ICM is Pic(R) cup Pic(O)
    assert forall{i:i,j in [1..#classes]|(i eq j) eq IsIsomorphic(classes[i],classes[j] : UseSpecializedMethod:=true)};

    "-------------------------------------------------------------";
    "-------------------------------------------------------------";
    "A quick test, with non trivial Class group, first with V=K and then with V = K^2";
    // Pic(O) is non trivial
    AttachSpec(path_to_package cat "AlgEt/spec");
    AttachSpec(path_to_package cat "AlgEt/specMod");
    AttachSpec(path_to_package cat "AlgEt/specMtrx");
    SetVerbose("IsomModules",1);
    _<x>:=PolynomialRing(Integers()); 
    g:=x^2+15;
    K:=EtaleAlgebra(g);
    nf:=Components(K);
    pi:=PrimitiveElement(K);
    R:=Order([pi]);
    m:=NaturalAction(K,K);
    time icm0:=ICM(R);
    time icm:=IsomorphismClasses(R,m : Method:="Magma", UseSpecializedMethod:=true);
    assert #icm eq #icm0;
    assert forall{i:i,j in [1..#classes]|(i eq j) eq IsIsomorphic(classes[i],classes[j] : UseSpecializedMethod:=true)};
    time classes:=IsomorphismClasses(R,m : Method:="Magma");
    assert #classes eq #icm; // Since R is Bass, the ICM is Pic(R) cup Pic(O)
    assert forall{i:i,j in [1..#classes]|(i eq j) eq IsIsomorphic(classes[i],classes[j] : UseSpecializedMethod:=true)};
    time classes:=IsomorphismClasses(R,m : Method:="julia -J /tmp/Hecke.so " cat path_to_package cat "AlgEt/AlgEtQMod/"); // this is the ICM
    assert #classes eq #icm; // Since R is Bass, the ICM is Pic(R) cup Pic(O)
    assert forall{i:i,j in [1..#classes]|(i eq j) eq IsIsomorphic(classes[i],classes[j] : UseSpecializedMethod:=true)};
    // V = K^2
    V:=EtaleAlgebra(&cat[nf : i in [1..2]]);
    m:=NaturalAction(K,V);
    time classes:=IsomorphismClasses(R,m : Method:="Magma");
    assert #classes eq 6; // because R is Bass 
    time classes:=IsomorphismClasses(R,m : Method:="julia -J /tmp/Hecke.so " cat path_to_package cat "AlgEt/AlgEtQMod/");
    assert #classes eq 6; // because R is Bass 
    assert forall{i:i,j in [1..#classes]|(i eq j) eq IsIsomorphic(classes[i],classes[j] : UseSpecializedMethod:=true)};

    "-------------------------------------------------------------";
    "-------------------------------------------------------------";
    "A very quick test, where V = K, that is we compute the ICM";
    AttachSpec(path_to_package cat "AlgEt/spec");
    AttachSpec(path_to_package cat "AlgEt/specMod");
    AttachSpec(path_to_package cat "AlgEt/specMtrx");
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
    time classes:=IsomorphismClasses(R,m : Method:="julia -J /tmp/Hecke.so " cat path_to_package cat "AlgEt/AlgEtQMod/", UseSpecializedMethod:=false);
    time icm:=ICM(R);
    assert #classes eq #icm; // Since R is Bass, the ICM is Pic(R) cup Pic(O)
    assert forall{i:i,j in [1..#classes]|(i eq j) eq IsIsomorphic(classes[i],classes[j] : UseSpecializedMethod:=true)};
    time classes:=IsomorphismClasses(R,m); // using icm
    assert #classes eq #icm; // Since R is Bass, the ICM is Pic(R) cup Pic(O)
    assert forall{i:i,j in [1..#classes]|(i eq j) eq IsIsomorphic(classes[i],classes[j] : UseSpecializedMethod:=true)};

    "-------------------------------------------------------------";
    "-------------------------------------------------------------";
    "A test that takes ~3 minutes, where V = K1^2 x K2";
    AttachSpec(path_to_package cat "AlgEt/spec");
    AttachSpec(path_to_package cat "AlgEt/specMod");
    AttachSpec(path_to_package cat "AlgEt/specMtrx");
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
    time classes:=IsomorphismClasses(R,m : Method:="julia -J /tmp/Hecke.so " cat path_to_package cat "AlgEt/AlgEtQMod/"); // changes this line accordingly to wheter you have used Hecke.Build() or not,
                                                                                                       // and to the appriopriate path to the the packages AlgEtQ
    assert #classes eq 4;

    // ########## 
    // Tests specific for the Power of Bass cases
    // #########

    "-------------------------------------------------------------";
    "-------------------------------------------------------------";
    "--------------------------PowerBass.m------------------------";
    "-------------------------------------------------------------";
    "-------------------------------------------------------------";

    AttachSpec("~/AlgEt/spec");
    AttachSpec("~/AlgEt/specMod");
    AttachSpec("~/AlgEt/specMtrx");
    
    P<x>:=PolynomialRing(Integers());
    f:=x^2 + x + 2;
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    R:=Order([2*F]);

    O:=MaximalOrder(A);
    ff:=Conductor(R);
    V:=DirectProduct([A,A,A]);
    m:=NaturalAction(A,V);
    Vnf:=Components(V);
    Anfpoly:=[ DefiningPolynomial(L) : L in Components(A) ];
    Vnfpoly:=[ DefiningPolynomial(L) : L in Vnf ];
    MO:=Module(R,m,<1*MaximalOrder(Vnf[i]) : i in [1..#Vnf]>);
    ff:=O!!Conductor(R);
    test,ff_prod:=IsProductOfIdeals(ff);
    assert test;
    ind:=[ #Vnfpoly+1 - Index(Reverse(Vnfpoly),fA)  : fA in Anfpoly ]; // last occurence of each 
                                                                       // number field from the dec of K 
                                                                       // in the decomposition of V
    Mff:=Module(R,m,<ff_prod[Index(Anfpoly,Vnfpoly[i])] : i in [1..#Vnf]>);
    candidates:=MaximalIntermediateModules(MO,Mff);
    for M in candidates do
        _:=DirectSumRepresentation(M);
        _:=SteinitzClass(M);
    end for;

    ///////////////////

    AttachSpec("~/AlgEt/spec");
    AttachSpec("~/AlgEt/specMod");
    AttachSpec("~/AlgEt/specMtrx");

    // takes a while
    P<x>:=PolynomialRing(Integers());
    f:=x^4 + 3*x^3 + 8*x^2 + 39*x + 169;
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    q:=Integers() ! (Coefficients(f)[1]^(2/Degree(f)));
    R:=Order([F,q/F]);
    ver:=[62,97,144,206,286];
    res:=[ #IsomorphismClassesOverBassOrder(R,i) : i in [1..5] ];
    assert res eq ver;
    
    ///////////////////
    
    AttachSpec("~/AlgEt/spec");
    AttachSpec("~/AlgEt/specMod");
    AttachSpec("~/AlgEt/specMtrx");

    P<x>:=PolynomialRing(Integers());
    f:=x^6-x^5+2*x^4-2*x^3+4*x^2-4*x+8;
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    R:=Order([F,ComplexConjugate(F)]);
    iso_cl:=IsomorphismClassesOverBassOrder(R,3);
    assert #iso_cl eq 6;

    M:=iso_cl[1];
    
    _:=[ StandardDirectSumRepresentation(M) : M in iso_cl ];

    t,s:=IsIsomorphicOverBassOrder(iso_cl[1],iso_cl[1]);
    assert t;
    t,s:=IsIsomorphicOverBassOrder(iso_cl[1],iso_cl[2]);
    assert not t;

    for M1,M2 in iso_cl do
        t,s:=IsIsomorphicOverBassOrder(M1,M2);
        assert t eq IsIsomorphic(M1,M2 : UseSpecializedMethod:=false ); // generic test
    end for;

    // ########## 
    // Slow Tests
    // #########

    "-------------------------------------------------------------";
    "-------------------------------------------------------------";
    "------------------IsomModules.m Slow Tests-------------------";
    "-------------------------------------------------------------";
    "-------------------------------------------------------------";
    "The following test should require around ~10000 for the julia sieving. Here R is Bass, but we compute it with the slow method to test it.";
    // Pic(O) is non triviali. Class construction is rther fast, isomorphism sieveng requires ~10000 secs using julia.
    AttachSpec(path_to_package cat "AlgEt/spec");
    AttachSpec(path_to_package cat "AlgEt/specMod");
    AttachSpec(path_to_package cat "AlgEt/specMtrx");
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
    time classes:=IsomorphismClasses(R,m : Method:="julia -J /tmp/Hecke.so " cat path_to_package cat "AlgEt/AlgEtQMod/");
    assert #classes eq 8; // because R is Bass 

    AttachSpec(path_to_package cat "AlgEt/spec");
    AttachSpec(path_to_package cat "AlgEt/specMod");
    AttachSpec(path_to_package cat "AlgEt/specMtrx");
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
    time classes:=IsomorphismClasses(R,m : Method:="julia -J /tmp/Hecke.so " cat path_to_package cat "AlgEt/AlgEtQMod/");
    assert #classes eq 6; // Example 6.1 in "Computing abelian varieties over finite fields isogenous to a power", by Marseglia


    
    "The following test is much bigger! (so far it never finished)";
    // a much bigger test!
    AttachSpec(path_to_package cat "AlgEt/spec");
    AttachSpec(path_to_package cat "AlgEt/specMod");
    AttachSpec(path_to_package cat "AlgEt/specMtrx");
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
    " The situation is the following: (we denote inclusion by c)
    ZFV c S c T c O,
    each inclusion with index 2.
    S has CohenMacaulayType 2, while ZFV,T and O are Gorenstein.
    All have only one prime above 2.
    All have trivial PicardGroup
    If ZFV were Bass then we would have the following number of isomorphism classes of ZFV modules in V=K^2:";
    ss:={1,2,3,4};
    #[ c : c in car<ss,ss,ss,ss> | c[1] le c[2] and c[2] le c[3] and c[3] le c[4] ];
    // we would get 35 classes.
    "But ZFV is NOT Bass";
    assert #PicardGroup(ZFV) eq 1;

    O:=MaximalOrder(K);
    nf:=Components(K);
    V:=EtaleAlgebra(&cat[nf : i in [1..2]]);
    m:=NaturalAction(K,V);
    time classes:=IsomorphismClasses(ZFV,m : Method:="julia -J /tmp/Hecke.so " cat path_to_package cat "AlgEt/AlgEtQMod/");
    "We get this many classes:";
    #classes;
   



















