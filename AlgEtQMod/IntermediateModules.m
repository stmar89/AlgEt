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

freeze;

declare verbose IntermediateModules, 2;

intrinsic MinimalIntermediateModules(I::AlgEtQMod,J::AlgEtQMod)->SetIndx[AlgEtQMod]
{Given S-modules J subset I, returns the minimal (with respect to inclusion) S-modules M such that J subset M subset I.}
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
        for p in PrimeDivisors(#Q) do
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

intrinsic IntermediateModules(I::AlgEtQMod,J::AlgEtQMod)->SetIndx[AlgEtQMod]
{Given fractional S-ideals J subset I, returns all the fractional S-ideals K such that J subset K subset I. They are produced recursively using from the minimal ones.}
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

intrinsic MaximalIntermediateModules(I::AlgEtQMod,J::AlgEtQMod)->SetIndx[AlgEtQMod]
{Given S-modules J subset I, returns the maximal (with respect to inclusion) S-modules K such that J subset K subset I.}
    assert2 J subset I; // "the ideal J needs to be inside I";
    S:=Order(I);
    assert2 S eq Order(J); // "The ideals must be over the same order";
    if J eq I then 
        return {@ J @};
    else
        // If I c K c I is maximal then I/K is a simple module, that is, isomorphic to S/P for some prime P of S.
        // In particular, if p is the characteristic of S/P, then p*I c K c I.
        // Note that I/pI is an Fp-vector space.
        _,m:=UniverseAlgebra(I);
        gens_J_over_S:=Generators(J);
        Q,q:=Quotient(I,J);
        max_submod:={@ @};
        for p in PrimeDivisors(#Q) do
            mp:=[ p*(Q.j) : j in [1..Ngens(Q)] ];
            mp_in_I:=[ x@@q : x in mp  ] cat gens_J_over_S;
            Qp,qp:=quo<Q|mp>;
            rk:=Ngens(Qp);
            assert #Qp eq p^rk;
            Fp:=FiniteField(p);
            // Note: Qp = Fp^rk as a vector space
            // So its S-module structure can be realized over Fp.
            // More precisely, given generators ai of S over Z, 
            // to find the (minimal) sub-S-modules, it is enough to look at Qp as module over Fp[ai]
            // In matrices we store the matrix representations of the action of the ZBasis of S on Qp
            matrices:=[Matrix(Fp,[ Eltseq(qp(q(m(zS)*(Qp.j@@qp@@q)))) : j in [1..Ngens(Qp)]]) : zS in ZBasis(S) ];
            matrices:=Setseq(Seqset(matrices)); //possibly there are repetiions.
            Qp_Rmod:=RModule(matrices);
            max_Rmod:=MaximalSubmodules(Qp_Rmod);
            max_submod join:={@ Module(S,m, mp_in_I cat 
                                    [(Qp!Eltseq(Qp_Rmod!b))@@qp@@q : b in Basis(max)]) : max in max_Rmod @};
        end for;
        return max_submod;
    end if;
end intrinsic;

intrinsic IntermediateModulesWithTrivialExtension(I::AlgEtQMod,J::AlgEtQMod,O::AlgEtQOrd)->SetIndx[AlgEtQMod]
{Given S-modules J subset I, and overorder O of S, it returns all the S-modules N such that J subset N subset I and NO=I. Note: we need O!!I eq I. They are produced recursively using from the maximal ones}
    V,m:=UniverseAlgebra(I);
    VJ,mJ:=UniverseAlgebra(J);
    S:=Order(I);
    require S eq Order(J) : "The modules must be over the same order";
    require V eq VJ and forall{b : b in Basis(Algebra(S)) | m(b) eq mJ(b)} : "the modules are not compatible";
    require J subset I : "The module J needs to be inside I";
    require S subset O : "The order O must be bigger than the order of I"; 
    IO:=O!!I;
    require S!!IO eq I : "I must be an O-module";
    queue:={@ I @};
    output:={@ I @};
    done:={@ @};
    while #queue gt 0 do
        pot_new:=&join[MaximalIntermediateModules(elt,J) : elt in queue ];
        pot_new:={@ K : K in pot_new | O!!K eq IO @}; // only the ones with trivial extension.
        output join:=pot_new;
        done join:=queue;
        // Note: if O!!K is not IO then every submodule of K will also not have trivial extension. Hence such K's can be excluded from the recursion.
        queue := pot_new diff done;
    end while;
    return output;
end intrinsic;

/* TESTS
    // see all_tests_AlgEtQMod.m

*/
