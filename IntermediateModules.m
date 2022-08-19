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
    require Order(I) eq Order(J) : "The ideals must be over the same order";
    require V eq VJ and forall{b : b in Basis(Algebra(Order(I))) | m(b) eq mJ(b)} : "the modules are not compatible";
    require J subset I : "The ideal J needs to be inside I";
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
    time #MinimalIntermediateModules(NE,N2);
    time #IntermediateModules(NE,N2);

*/
