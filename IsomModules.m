/* vim: set syntax=magma :*/

freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose IsomModules, 2;

/* TODO
    - add how to use julia to the description
    - TESTS: squarefree with non trivial class group
             power of Bass non trivial class group
    - when looping over the different k's, do the classes_k's all contain the same number of objects?
      see CONJECTURE BELOW
*/

//------------
// Isomorphism Testing for modules
//------------

intrinsic IsIsomorphic(I::AlgEtMod,J::AlgEtMod : Method:="julia") -> BoolElt
{Given two modules I and J returns wheater they are isomorphic.
The vararg Method allows to choose if the isomorphism testing is done with "Magma", very slow, or with the Hecke/Nemo package for julia, which is much faster.
In the latter case, the Method should be of the form
"julia path/to/AlgEt/", or if Hecke/Nemo has been built,"julia -J /tmp/Hecke.so ~/path/to/AlgEt/" (the ".so" might be different according to your SO. See the documentation of Hecke.Build()).}
    require Method eq "Magma" or Method[1..5] eq "julia" : "Method should be Magma or julia or a string of julia command";
    V,m:=UniverseAlgebra(I);
    VJ,mJ:=UniverseAlgebra(I);
    S:=Order(I);
    require S eq Order(J) : "The ideals must be over the same order";
    require V eq VJ and forall{b : b in Basis(Algebra(S)) | m(b) eq mJ(b)} : "the modules are not compatible";
   
    //require EquationOrder(Algebra(S)) subset S : "Implemented only for modules over orders containing the equation order";
    pi:=PrimitiveElement(Algebra(S));
    if not pi in S then
        pi:=pi*Exponent(Quotient(pi*S + OneIdeal(S),OneIdeal(S)));
    end if;
    pi:=m(pi);
    // pi is now a primitive element of K such that Z[pi] c S.
    
    // based on the following, two Z[pi] modules are isomorphic iff the matrices representing multiplcition by pi are Z-conjugte 
    matI:=ChangeRing(Matrix(AbsoluteCoordinates([pi*z : z in ZBasis(I)],ZBasis(I))),Integers());
    matJ:=ChangeRing(Matrix(AbsoluteCoordinates([pi*z : z in ZBasis(J)],ZBasis(J))),Integers());
    if Method eq "Magma" then
        test:=AreGLConjugate(matI,matJ);
    else //Method eq "julia ...."
        ID:=&cat[ Sprint(Random(0,9)) : i in [1..20] ];
        file:="tmp_" cat ID cat ".txt";
        str:=Sprintf("[\n%o,\n%o\n]\n", Eltseq(matI),Eltseq(matJ));
        fprintf file, "%o",str;
        indices:=eval(Pipe(Method cat "IsIsomorphic_julia_script.jl " cat file,""));
        if #indices eq 2 then
        // the matrices are not conjugate
            test:=false;
        else
            assert #indices eq 1;
            test:=true;
        end if;
        Pipe("rm " cat file cat " || true ",""); 
    end if;
    return test;
end intrinsic;

//------------
// Isomorphism Classes for modules
//------------

intrinsic IsomorphismClasses(R::AlgEtOrd,m::Map : Method:="julia") -> SeqEnum[AlgEtMod]
{Given an order R in some AlgEt K, where K acts on some V, by m:K->V, returns representatives of hte isomorphism classes of the S-module lattices in V.
The vararg Method allows to choose if the isomorphism testing is done with "Magma", very slow, or with the Hecke/Nemo package for julia, which is much faster.
In the latter case, the Method should be of the form "julia path/to/AlgEt/", or if Hecke/Nemo has been built,"julia -J /tmp/Hecke.so ~/path/to/AlgEt/" (the ".so" might be different according to your SO. See the documentation of Hecke.Build()).}
    require Method eq "Magma" or Method[1..5] eq "julia" : "Method should be Magma or julia or a string of julia command";
    V:=Codomain(m);
    K:=Domain(m);
    pi:=PrimitiveElement(K);
    if not pi in R then
        pi:=pi*Exponent(Quotient(pi*R + OneIdeal(R),OneIdeal(R)));
    end if;
    pi:=m(pi);
    O:=MaximalOrder(K);
    ff:=Conductor(R);
    Vnf:=NumberFields(V);
    Knfpoly:=[ DefiningPolynomial(L) : L in NumberFields(K) ];
    Vnfpoly:=[ DefiningPolynomial(L) : L in Vnf ];
    MO:=Module(R,m,<1*MaximalOrder(Vnf[i]) : i in [1..#Vnf]>);
    ff:=O!!Conductor(R);
    test,ff_prod:=IsProductOfIdeals(ff);
    assert test;
    ind:=[ #Vnfpoly+1 - Index(Reverse(Vnfpoly),fK)  : fK in Knfpoly ]; // last occurence of each number field from the dec of K 
                                                                       // in the decomposition of V
    Mff:=Module(R,m,<ff_prod[Index(Knfpoly,Vnfpoly[i])] : i in [1..#Vnf]>);
    vprint IsomModules,1 : "candidates:";
    vtime IsomModules,1 : candidates:=IntermediateModulesWithTrivialExtension(MO,Mff,O);
    vprint IsomModules,1 : #candidates;
    PO,pO:=PicardGroup(O);
    classes:=[];
    for g in PO do
        if g eq Zero(PO) then
            candidates_k:=candidates; 
        else
            _,Ik:=CoprimeRepresentative(pO(g),ff); // Ik+ff=O, Ik cap ff = Ik*ff
            test,Ik_prod:=IsProductOfIdeals(Ik);
            assert test;
            MffIk:=Module(R,m,<i in ind select 
                                            ff_prod[Index(Knfpoly,Vnfpoly[i])] 
                                        else ff_prod[Index(Knfpoly,Vnfpoly[i])]*Ik_prod[Index(Knfpoly,Vnfpoly[i])] 
                                                : i in [1..#Vnf]>);
            e:=ChineseRemainderTheorem(ff,Ik,One(K),Zero(K)); // e in Ik, e-1 in ff.
            assert not IsZeroDivisor(e);
            e:=Components(e);
            e:=V!< i in ind select e[Index(Knfpoly,Vnfpoly[i])] else Vnf[i] ! 1 : i in [1..#Vnf] >;
            ik:=map<V->V | x:->x*e >;
            // ik induces the isomorphism between (O1^s1 + ... + On^sn)/(f1^s1 + ... + fn^sn),
            // and (O1^(s1-1)+I1 + ... + On^(sn-1)+In)/(f1^(s1-1)+f1I1 + ... + fn^(sn-1)+fnIn),
            // where Ik = I1 + .... + In.
            candidates_k:=[ Module(R,m, [ik(z) : z in ZBasis(M)]) : M in candidates];
        end if;
        if Method eq "Magma" then
            classes_k:=[];
            for M in candidates_k do
                if not exists{ N : N in classes_k | IsIsomorphic(M,N : Method:=Method) } then
                    Append(~classes_k,M);
                end if;
            end for;
        else //Method eq "julia ...."
            ID:=&cat[ Sprint(Random(0,9)) : i in [1..20] ];
            file:="tmp_" cat ID cat ".txt";
            str:="[\n";
            for i->I in candidates do
                mat:=Matrix(AbsoluteCoordinates([pi*z : z in ZBasis(I)],ZBasis(I)));
                str cat:=Sprintf("%o,\n", Eltseq(mat));
            end for;
            str:=Prune(Prune(str)) cat "\n]\n";
            fprintf file, "%o",str;
            indices_of_classes_k:=eval(Pipe(Method cat "IsIsomorphic_julia_script.jl " cat file,""));
            vprintf IsomModules,1 : "indices_of_classes_k = %o\n",indices_of_classes_k;
            classes_k:=[candidates_k[i] : i in indices_of_classes_k];
            Pipe("rm " cat file cat " || true ",""); 
        end if;
        Append(~classes,classes_k);
    end for;
    assert forall{ c : c in classes | #c eq #classes[1] }; // CONJECTURE THEY ALL CONTAIN THE SAME NUMBER OF CLASSES
    classes:=&cat(classes);
    return classes;
end intrinsic;

/* TEST
  
    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/Modules.m");
    Attach("~/packages_github/AlgEt/IntermediateModules.m");
    Attach("~/packages_github/AlgEt/IsomModules.m");
    SetVerbose("IsomModules",1);
    _<x>:=PolynomialRing(Integers()); 
    g:=x^2+5;
    K:=EtaleAlgebra(g);
    nf:=NumberFields(K);
    pi:=PrimitiveElement(K);
    R:=Order([7*pi]);
    m:=NaturalAction(K,K);
    time classes:=IsomorphismClasses(R,m : Method:="julia -J /tmp/Hecke.so ~/packages_github/AlgEt/"); // this is the ICM
    time icm:=ICM(R);
    assert #classes eq #icm; // Since R is Bass, the ICM is Pic(R) cup Pic(O)

    // Pic(O) is non trivial
    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/Modules.m");
    Attach("~/packages_github/AlgEt/IntermediateModules.m");
    Attach("~/packages_github/AlgEt/IsomModules.m");
    SetVerbose("IsomModules",1);
    _<x>:=PolynomialRing(Integers()); 
    g:=x^2+15;
    K:=EtaleAlgebra(g);
    nf:=NumberFields(K);
    pi:=PrimitiveElement(K);
    R:=Order([pi]);
    m:=NaturalAction(K,K);
    time classes:=IsomorphismClasses(R,m : Method:="julia -J /tmp/Hecke.so ~/packages_github/AlgEt/"); // this is the ICM
    time icm:=ICM(R);
    assert #classes eq #icm; // Since R is Bass, the ICM is Pic(R) cup Pic(O)
    // V = K^2
    V:=EtaleAlgebra(&cat[nf : i in [1..2]]);
    m:=NaturalAction(K,V);
    time classes:=IsomorphismClasses(R,m : Method:="julia -J /tmp/Hecke.so ~/packages_github/AlgEt/");
    assert #classes eq 6; // because R is Bass 
    // THERE IS A BUG in the julia code :  V = K^3
    V:=EtaleAlgebra(&cat[nf : i in [1..3]]);
    m:=NaturalAction(K,V);
    time classes:=IsomorphismClasses(R,m : Method:="julia -J /tmp/Hecke.so ~/packages_github/AlgEt/");
    assert #classes eq 8; // because R is Bass 

    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/Modules.m");
    Attach("~/packages_github/AlgEt/IntermediateModules.m");
    Attach("~/packages_github/AlgEt/IsomModules.m");
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
    time classes:=IsomorphismClasses(R,m : Method:="julia -J /tmp/Hecke.so ~/packages_github/AlgEt/"); // changes this line accordingly to wheter you have used Hecke.Build() or not,
                                                                                                       // and to the appriopriate path to the the packages AlgEt
    assert #classes eq 4;

    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/Modules.m");
    Attach("~/packages_github/AlgEt/IntermediateModules.m");
    Attach("~/packages_github/AlgEt/IsomModules.m");
    SetVerbose("IsomModules",1);
    _<x>:=PolynomialRing(Integers()); 
    g:=x^6-x^5+2*x^4-2*x^3+4*x^2-4*x+8;
    K:=EtaleAlgebra(g);
    nf:=NumberFields(K);
    q:=2;
    pi:=PrimitiveElement(K);
    R:=Order([pi,q/pi]);
    O:=MaximalOrder(K);
    assert IsBass(R);
    m:=NaturalAction(K,K);
    time classes:=IsomorphismClasses(R,m : Method:="julia -J /tmp/Hecke.so ~/packages_github/AlgEt/"); // this is the ICM
    time icm:=ICM(R);
    assert #classes eq #icm; // Since R is Bass, the ICM is Pic(R) cup Pic(O)

    // THIS LAST TEST IS VERY SLOW. Computing the candidates (1774) takes 470 secs. The julia routine to sieve out the representatives does not finish in one night.
    // V = K^3
    V:=EtaleAlgebra(&cat[nf : i in [1..3]]);
    m:=NaturalAction(K,V);
    time classes:=IsomorphismClasses(R,m : Method:="julia -J /tmp/Hecke.so ~/packages_github/AlgEt/");
    assert #classes eq 6; // Example 6.1 in "Computing abelian varieties over finite fields isogenous to a power", by Marseglia

*/
