/////////////////////////////////////////////////////
// Stefano Marseglia, stefano.marseglia89@gmail.com
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

declare verbose IntermediateIdeals, 2;

///# Intermediate fractional ideals
/// Let $I$ and $J$ be orders or fractional ideals in an étale algebra over $\mathbb{Q}$ satisfying $J subseteq I$.
/// Since the quotient $I/J$ is finite, there are finitely many fractional ideals $K$ (over a fixed order) such that $J \subseteq K \subseteq I$.
/// The following intrinsic allow to compute them.

/// Given fractional $S$-ideals $I$ and $J$ such that $J \subseteq I$, returns the minimal (with respect to inclusion) fractional $S$-ideals $K$ such that $J \subsetneq K \subseteq I$. Note that $J$ is never in the output.
intrinsic MinimalIntermediateIdeals(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]
{Given fractional S-ideals I and J such that J subseteq I, returns the minimal (with respect to inclusion) fractional S-ideals K such that J subsetneq K subseteq I. Note that J is never in the output.}
    assert2 J subset I; // "the ideal J needs to be inside I";
    S:=Order(I);
    assert2 S eq Order(J); // "The ideals must be over the same order";
    if J eq I then 
        return {@ @};
    else
        gens_J_over_S:=Generators(J);
        Q,q:=Quotient(I,J);
        min_ideals:={@ @};
        // Following code is based on the following Lemma.
        // LEMMA: Let M be a finite non-zero S-module (eg I/J). Let N be a minimal submodule of M (non-zero, possibly N=M).
        // Then there exists a prime number p dividing #M such that pN=0.
        // In particular, N is a submodule of the kernel Mp of the multiplication-by-p map on M.
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
            matrices:=[Matrix(Fp,[ Eltseq(Qp!(q(zS*(Qp.j@@q)))) : j in [1..Ngens(Qp)]]) : zS in ZBasis(S) ];
            matrices:=Setseq(Seqset(matrices)); //possibly there are repetiions.
            Qp_Rmod:=RModule(matrices);
            min_Rmod:=MinimalSubmodules(Qp_Rmod);
            min_ideals join:={@ Ideal(S, gens_J_over_S cat 
                                    [(Q!(Qp!Eltseq(Qp_Rmod!b)))@@q : b in Basis(min)]) : min in min_Rmod @};
        end for;
        return min_ideals;
    end if;
end intrinsic;

/// Given fractional $S$-ideals $I$ and $J$ such that $J$ \subseteq $I$, returns all the fractional $S$-ideals $K$ such that $J \subseteq K \subseteq I$. They are produced recursively from the minimal ones. The output includes $I$ and $J$.
intrinsic IntermediateIdeals(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]
{Given fractional S-ideals I and J such that J subseteq I, returns all the fractional S-ideals K such that J subseteq K subseteq I. They are produced recursively from the minimal ones. The output includes I and J.}
    require J subset I : "The ideal J needs to be inside I";
    require Order(I) eq Order(J) : "The ideals must be over the same order";
    queue:={@ J @};
    output:={@ J @};
    done:={@ @};
    while #queue gt 0 do
        pot_new:=&join[MinimalIntermediateIdeals(I,elt) : elt in queue ];
        output join:=pot_new;
        done join:=queue;
        queue := pot_new diff done;
    end while;
    return output;
end intrinsic;

/// Given fractional $S$-ideals $I$ and $J$ such that $J \subseteq I$, returns all the fractional $S$-ideals $K$ such that $J \subseteq K \subseteq I$ and $(K:K)=S$. They are produced recursively from the minimal ones. The output includes $I$, if $(I:I)=S$, and $J$, if $(J:J)=S$.
intrinsic IntermediateIdealsWithPrescribedMultiplicatorRing(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]
{Given fractional S-ideals I and J such that J subseteq I, returns all the fractional S-ideals K such that J subseteq K subseteq I and (K:K)=S. They are produced recursively from the minimal ones. The output includes I, if (I:I)=S, and J, if (J:J)=S.}
    require J subset I : "The ideal J needs to be inside I";
    S:=Order(I);
    require S eq Order(J) : "The ideals must be over the same order";
    queue:={@ J @};
    if MultiplicatorRing(J) eq S then
        output:={@ J @};
    else
        output:={@  @};
    end if;
    done:={@ @};
    while #queue gt 0 do
        pot_new:=&join[MinimalIntermediateIdeals(I,elt) : elt in queue ];
        output join:={@ K : K in pot_new | not K in done and MultiplicatorRing(K) eq S @};
        done join:=queue;
        queue := pot_new diff done;
    end while;
    return output;
end intrinsic;

/// Given fractional $S$-ideals $I$ and $J$ such that $J \subseteq I$, returns the maximal (with respect to inclusion) fractional $S$-ideals $K$ such that $J \subseteq K \subsetneq I$. Note that $I$ is never in the output, while $J$ is in the output if and only if the $S$-module $I/J$ is simple, in which case the output consists only of $J$.
intrinsic MaximalIntermediateIdeals(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]
{Given fractional S-ideals I and J such that J subseteq I, returns the maximal (with respect to inclusion) fractional S-ideals K such that J subseteq K subsetneq I. Note I is never in the output, while J is in the output if and only if the S-module I/J is simple, in which case the output consists only of J.}
    assert2 J subset I; // "the ideal J needs to be inside I";
    S:=Order(I);
    assert2 S eq Order(J); // "The ideals must be over the same order";
    if J eq I then 
        return {@ @};
    else
        // If J c K c I is maximal then I/K is a simple module, that is, isomorphic to S/P for some prime P of S.
        // In particular, if p is the characteristic of S/P, then p*I c K c I.
        // Note that I/pI is an Fp-vector space.
        Q,q:=Quotient(I,J);
        gens_J_over_S:=Generators(J);
        max_ideals:={@ @};
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
            matrices:=[Matrix(Fp,[ Eltseq(qp(q(zS*(Qp.j@@qp@@q)))) : j in [1..Ngens(Qp)]]) : zS in ZBasis(S) ];
            matrices:=Setseq(Seqset(matrices)); //possibly there are repetiions.
            Qp_Rmod:=RModule(matrices);
            max_Rmod:=MaximalSubmodules(Qp_Rmod);
            max_ideals join:={@ Ideal(S, mp_in_I cat 
                                    [(Qp!Eltseq(Qp_Rmod!b))@@qp@@q : b in Basis(max)]) : max in max_Rmod @};
        end for;
        return max_ideals;
    end if;
end intrinsic;

/// Given fractional $S$-ideals $I$ and $J$ and an order $O$ such that 
/// - $S \subseteq O$,  
/// - $J \subseteq I$, and 
/// - $O \subseteq (I:I)$,
/// returns all the fractional $S$-ideals $K$ such that 
/// - $J \subseteq K \subseteq I$, and 
/// - $O\cdot K = I$. 
/// Note that the output always contains $I$. 
/// The output is produced by recursively computing maximal intermediate ideals.
intrinsic IntermediateIdealsWithTrivialExtension(I::AlgEtQIdl,J::AlgEtQIdl, O::AlgEtQOrd)->SetIndx[AlgEtQIdl]
{Given fractional S-ideals I and J and an order O such that 
- S subseteq O,  
- J subseteq I, and 
- O subseteq (I:I),
returns all the fractional S-ideals K such that 
- J subseteq K subseteq I, and 
- O!!K = I. 
Note that the output always contains I. The output is produced by recursively computing maximal intermediate ideals.}
    require J subset I : "The ideal J needs to be inside I";
    S:=Order(I);
    require S eq Order(J) : "The ideals must be over the same order";
    require S subset O : "O is not an overorder of Order(I)";
    IO:=O!!I;
    require S!!IO eq I : "I is not an O-ideal";
    queue:={@ I @};
    output:={@ I @};
    done:={@ @};
    while #queue gt 0 do
        pot_new:=&join[MaximalIntermediateIdeals(elt,J) : elt in queue ];
        pot_new:={@ K : K in pot_new | O!!K eq IO @}; //we keep only the ones with trivial extension
        output join:={@ K : K in pot_new | not K in done @};
        done join:=queue;
        // Note: if O!!K is not IO, then all the submodules of K will also not have trivial extension IO.
        // Hence we don't need to continue the recursion on K.
        queue := pot_new diff done; 
    end while;
    return output;
end intrinsic;

/// Given fractional $S$-ideals $I$ and $J$ and an order $O$ such that 
/// - $S \subseteq O$,  
/// - $J \subseteq I$, and 
/// - $O \subseteq (I:I)$, 
/// returns all the fractional $S$-ideals $K$ satisfying
/// - $J \subseteq K \subseteq I$, 
/// - $O\cdot K = I$, and 
/// - $(K:K) = S$. 
/// In particular, the output contains $I$ if and only if $O = (I:I) = S$. The output is produced by recursively computing maximal intermediate ideals.
intrinsic IntermediateIdealsWithTrivialExtensionAndPrescribedMultiplicatorRing(I::AlgEtQIdl,J::AlgEtQIdl, O::AlgEtQOrd)->SetIndx[AlgEtQIdl]
{Given fractional S-ideals I and J and an order O such that 
- S subseteq O,  
- J subseteq I, and 
- O subseteq (I:I), 
returns all the fractional S-ideals K satisfying
- J subseteq K subseteq I, 
- O!!K = I, and 
- (K:K) = S. 
In particular, the output contains I if and only if O = (I:I) = S. The output is produced by recursively computing maximal intermediate ideals.}
    require J subset I : "The ideal J needs to be inside I";
    S:=Order(I);
    require S eq Order(J) : "The ideals must be over the same order";
    require S subset O : "O is not an overorder of Order(I)";
    IO:=O!!I;
    require MultiplicatorRing(I) subset O : "I is not an O-ideal"; 
    queue:={@ I @};
    if MultiplicatorRing(I) eq S then
        output:={@ I @};
    else
        output:={@ @};
    end if;
    done:={@ @};
    while #queue gt 0 do
        pot_new:=&join[MaximalIntermediateIdeals(elt,J) : elt in queue ];
        pot_new:={@ K : K in pot_new | O!!K eq IO @}; //we keep only the ones with trivial extension
        output join:={@ K : K in pot_new | not K in done and MultiplicatorRing(K) eq S @};
        done join:=queue;
        // Note: if O!!K is not IO, then all the submodules of K will also not have trivial extension IO.
        // Hence we don't need to continue the recursion on K.
        queue := pot_new diff done; 
    end while;
    return output;
end intrinsic;

/// Given fractional $S$-ideals $I$ and $J$ satisfying $J \subseteq I$, and a positive integer $N$, returns all the fractional $S$-ideals $K$ such that:
/// - $J \subseteq K \subseteq I$, and 
/// - $[I:K]=N$. 
/// The output is produced by recursively computing the maximal ones.
intrinsic IntermediateIdealsOfIndex(I::AlgEtQIdl,J::AlgEtQIdl,N::RngIntElt)->SetIndx[AlgEtQIdl]
{Given fractional S-ideals I and J satisfying J subseteq I, and a positive integer N, returns all the fractional S-ideals K such that 
- J subseteq K subseteq I, and 
- [I:K]=N. 
The output is produced by recursively computing the maximal ones.}
    require J subset I : "The ideal J needs to be inside I";
    require N gt 0 : "N must be a strictly positive integer";
    S:=Order(I);
    require S eq Order(J) : "The ideals must be over the same order";
    // early exits
    if N eq 1 then
        return {@ I @};
    end if;
    if Index(I,J) eq N then
        return {@ J @};
    elif (Index(I,J) mod N) ne 0 then
        output:={@ Universe({@ I @}) | @}; //empty set
        return output;
    end if;
    // we start the recursion
    queue:={@ I @};
    output:={@ @};
    done:={@ @};
    while #queue gt 0 do
        pot_new:=&join[MaximalIntermediateIdeals(elt,J) : elt in queue ];
        output join:={@ K : K in pot_new | Index(I,K) eq N @};
        done join:=queue;
        queue := {@ M : M in pot_new diff done | Index(I,M) lt N @}; // If [I:M] ge N then for all K c M we have 
                                                                     // that [I:K]>N. Hence we don't want such M's 
                                                                     // in the queue.
    end while;
    return output;
end intrinsic;

/* TESTS

    printf "### Testing IntermediateIdeals:";
	SetAssertions(2);
	_<x>:=PolynomialRing(Integers());
    f:=x^3-100*x^2-100*x-100;
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    O:=MaximalOrder(K);
    EO:=E!!OneIdeal(O);
    ff:=Conductor(E);
    mm:=MinimalIntermediateIdeals(EO,ff);
    printf ".";
    assert forall{I:I in mm | {@ I @} eq MinimalIntermediateIdeals(I,ff)};
    printf ".";
    nn:=IntermediateIdeals(EO,ff);
    printf ".";
    assert mm subset nn;
    assert EO in nn;
    assert ff in nn;
    printf ".";
    nn:=IntermediateIdealsWithPrescribedMultiplicatorRing(EO,ff);
    assert OneIdeal(E) in nn;
    printf ".";
    nn:=IntermediateIdealsWithPrescribedMultiplicatorRing(OneIdeal(O),O!!ff);
    assert OneIdeal(O) in nn;
    assert O!!ff in nn;
    printf ".";
    nn:=IntermediateIdealsWithPrescribedMultiplicatorRing(EO,OneIdeal(E));
    assert OneIdeal(E) in nn;
    printf ".";
    nn:=IntermediateIdealsWithPrescribedMultiplicatorRing(OneIdeal(E),ff);
    assert OneIdeal(E) in nn;
    printf ".";
    mm:=MaximalIntermediateIdeals(EO,ff);
    assert forall{I:I in mm | {@ I @} eq MaximalIntermediateIdeals(EO,I)};
    printf ".";
    nn:=IntermediateIdealsWithTrivialExtension(EO,ff,O);
    assert EO in nn;
    assert OneIdeal(E) in nn;
    printf ".";
    nn:=IntermediateIdealsWithTrivialExtensionAndPrescribedMultiplicatorRing(EO,ff,O);
    assert OneIdeal(E) in nn;
    printf ".";
    nn:=IntermediateIdealsWithTrivialExtensionAndPrescribedMultiplicatorRing(EO,OneIdeal(E),O);
    assert OneIdeal(E) in nn;
    printf ".";
    nn:=IntermediateIdealsWithTrivialExtensionAndPrescribedMultiplicatorRing(OneIdeal(O),O!!ff,O);
    assert OneIdeal(O) in nn;
    printf ".";
    f:=x^4+291*x^3-988*x^2-1000*x-1000;
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    O:=MaximalOrder(K);
    EO:=E!!OneIdeal(O);
    ff:=Conductor(E);
    mm:=MinimalIntermediateIdeals(EO,ff);
    printf ".";
    assert forall{I:I in mm | {@ I @} eq MinimalIntermediateIdeals(I,ff)};
    printf ".";
    nn:=IntermediateIdeals(EO,ff); //2 minutess
    printf ".";
    assert mm subset nn;
    assert EO in nn;
    assert ff in nn;
    printf ".";
    nn:=IntermediateIdealsWithPrescribedMultiplicatorRing(EO,ff);
    assert OneIdeal(E) in nn;
    printf ".";
    nn:=IntermediateIdealsWithPrescribedMultiplicatorRing(OneIdeal(O),O!!ff);
    assert OneIdeal(O) in nn;
    assert O!!ff in nn;
    printf ".";
    nn:=IntermediateIdealsWithPrescribedMultiplicatorRing(EO,OneIdeal(E));
    assert OneIdeal(E) in nn;
    printf ".";
    nn:=IntermediateIdealsWithPrescribedMultiplicatorRing(OneIdeal(E),ff);
    assert OneIdeal(E) in nn;
    printf ".";
    mm:=MaximalIntermediateIdeals(EO,ff);
    assert forall{I:I in mm | {@ I @} eq MaximalIntermediateIdeals(EO,I)};
    printf ".";
    nn:=IntermediateIdealsWithTrivialExtension(EO,ff,O);
    assert EO in nn;
    assert OneIdeal(E) in nn;
    printf ".";
    nn:=IntermediateIdealsWithTrivialExtensionAndPrescribedMultiplicatorRing(EO,ff,O);
    assert OneIdeal(E) in nn;
    printf ".";
    nn:=IntermediateIdealsWithTrivialExtensionAndPrescribedMultiplicatorRing(EO,OneIdeal(E),O);
    assert OneIdeal(E) in nn;
    printf ".";
    nn:=IntermediateIdealsWithTrivialExtensionAndPrescribedMultiplicatorRing(OneIdeal(O),O!!ff,O);
    assert OneIdeal(O) in nn;
    printf ".";
    SetAssertions(1);
    printf " all good!"; 

*/
