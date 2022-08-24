/* vim: set syntax=magma :*/

freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose IntermediateModules, 2;

/*TODO:

*/

intrinsic MinimalIntermediateModules(I::AlgEtMod,J::AlgEtMod)->SetIndx[AlgEtMod]
{ Given S-modules J subset I, returns the minimal (with respect to inclusion) S-modules M such that J subset M subset I. }
    assert2 J subset I; // "the ideal J needs to be inside I";
    S:=Order(I);
    assert2 S eq Order(J); // "The ideals must be over the same order";
    if J eq I then 
        return {@ @};
    else
        _,m:=UniverseAlgebra(I);
        gens_J_over_S:=Generators(J);
        Q,q:=Quotient(I,J);
        min_mods:={@ @};
        // Following code is based on the following Lemma.
        // LEMMA: Let M be a finite S-module (eg I/J). Let N be a minimal submodule of M. 
        // Then there exists a prime number p dividing #M such that pN=0.
        // In particular, N is a submoduel of the kernel Mp of the multiplication-by-p map on M.
        // PROOF: pN is a submodule of N. Since N is minimal pN=N or pN=0. 
        // The first case occurs if p does not divide the order of N.
        // The second occurs if p divides the order of N.
        // QED
        for p in PrimeDivisors(Index(I,J)) do
            mp:=hom<Q->Q | x:->p*x>;
            Qp:=Kernel(mp);
            rk:=Ngens(Qp);
            assert #Qp eq p^rk;
            Fp:=FiniteField(p);
            // Note: Qp = Fp^rk as a vector space
            // So its S-module structure can be realized over Fp.
            // More precisely, given generators ai of S over Z, 
            // to find the (minimal) sub-S-modules, it is enough to look at Qp as module over Fp[ai]
            // In matrices we store the matrix representations of the action of the ZBasis of S on Qp
            matrices:=[Matrix(Fp,[ Eltseq(Qp!(q(m(zS)*(Qp.j@@q)))) : j in [1..Ngens(Qp)]]) : zS in ZBasis(S) ];
            matrices:=Setseq(Seqset(matrices)); //possibly there are repetiions.
            Qp_Rmod:=RModule(matrices);
            min_Rmod:=MinimalSubmodules(Qp_Rmod);
            min_mods join:={@ Module(S,m, gens_J_over_S cat 
                                    [(Q!(Qp!Eltseq(Qp_Rmod!b)))@@q : b in Basis(min)]) : min in min_Rmod @};
        end for;
        return min_mods;
    end if;
end intrinsic;

intrinsic IntermediateModules(I::AlgEtMod,J::AlgEtMod)->SetIndx[AlgEtMod]
{ Given fractional S-ideals J subset I, returns all the fractional S-ideals K such that J subset K subset I. They are produced recursively using from the minimal ones }
    V,m:=UniverseAlgebra(I);
    VJ,mJ:=UniverseAlgebra(I);
    require Order(I) eq Order(J) : "The modules must be over the same order";
    require V eq VJ and forall{b : b in Basis(Algebra(Order(I))) | m(b) eq mJ(b)} : "the modules are not compatible";
    require J subset I : "The module J needs to be inside I";
    queue:={@ J @};
    output:={@ J @};
    done:={@ @};
    while #queue gt 0 do
        pot_new:=&join[MinimalIntermediateModules(I,elt) : elt in queue ];
        output join:=pot_new;
        done join:=queue;
        queue := pot_new diff done;
    end while;
    return output;
end intrinsic;

intrinsic IntermediateModulesWithPrescribedExtension(I::AlgEtMod,J::AlgEtMod,O::AlgEtOrd,M::AlgEtMod)->SetIndx[AlgEtMod]
{ Given S-modules J subset I, and overorder O of S and a fixed O-module M subset J, it returns all the S-modules N such that I subset N subset J and NO=M.
  They are produced recursively using from the minimal ones }
    V,m:=UniverseAlgebra(I);
    VJ,mJ:=UniverseAlgebra(I);
    S:=Order(I);
    require S eq Order(J) : "The modules must be over the same order";
    require V eq VJ and forall{b : b in Basis(Algebra(S)) | m(b) eq mJ(b)} : "the modules are not compatible";
    require J subset I : "The module J needs to be inside I";
    require S subset O : "The order O must be bigger than the order of I"; 
    require M subset O!!I : "The module M must be a subset of I";
    require Order(M) eq O : "The order of M must be O";
    queue:={@ J @};
    output:={@  @};
    done:={@ @};
    while #queue gt 0 do
        pot_new:=&join[MinimalIntermediateModules(I,elt) : elt in queue ];
        output join:={@ K : K in pot_new | not K in done and O!!K eq M @};
        done join:=queue;
        queue := pot_new diff done;
    end while;
    return output;
end intrinsic;

/* CONTINUE FROM HERE
intrinsic IntermediateIdealsWithPrescribedMultiplicatorRing(I::AlgEtIdl,J::AlgEtIdl)->SetIndx[AlgEtIdl]
{ Given fractional S-ideals J subset I, returns all the fractional S-ideals K such that (K:K)=S and  J subset K subset I. 
  They are produced recursively using from the minimal ones }
    require J subset I : "The ideal J needs to be inside I";
    S:=Order(I);
    require S eq Order(J) : "The ideals must be over the same order";
    queue:={@ J @};
    output:={@  @};
    done:={@ @};
    while #queue gt 0 do
        pot_new:=&join[MinimalIntermediateIdeals(I,elt) : elt in queue ];
        output join:={@ K : K in pot_new | not K in done and MultiplicatorRing(K) eq S @};
        done join:=queue;
        queue := pot_new diff done;
    end while;
    return output;
end intrinsic;
*/


/* TEST

    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/Modules.m");
    Attach("~/packages_github/AlgEt/IntermediateModules.m");
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

    Vnf,Ve:=NumberFields(V);
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
    time l:=IntermediateModules(NE,N2); #l; // 10 secs
    time l2:=IntermediateModulesWithPrescribedExtension(NE,N2,O,O!!NE); ~12 secs
    NEO:=O!!NE;
    assert #l2 eq #[ M : M in l | O!!M eq NEO ];

    // a much bigger test!

    AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/Modules.m");
    Attach("~/packages_github/AlgEt/IntermediateModules.m");
    Attach("~/packages_github/AlgEt/IsomModules.m");
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
    // But ZFV is NOT Bass!
    assert #PicardGroup(ZFV) eq 1;

    R:=ZFV;
    O:=MaximalOrder(K);
    nf:=NumberFields(K);
    V:=EtaleAlgebra(&cat[nf : i in [1..2]]);
    Vnf:=NumberFields(V);
    m:=NaturalAction(K,V);
    MO:=Module(R,m,<1*MaximalOrder(Vnf[i]) : i in [1..#Vnf]>);
    MOO:=O!!MO;
    ff:=O!!Conductor(R);
    test,ff_prod:=IsProductOfIdeals(ff);
    assert test;
    // the following works for this exmaple, not in general.
    Mff:=Module(R,m,<ff_prod[1] : i in [1..2]>);
    gensMff:=Generators(Mff);
    time l0:=IntermediateModulesWithPrescribedExtension(MO,Mff,O,MOO);
    time l1:=IntermediateModules(MO,Mff);
    time l2:=[ I : I in candidates | O!!I eq MO ];
    assert #l2 eq #l0;

*/
