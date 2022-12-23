/* vim: set syntax=magma :*/

freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose IntermediateIdeals, 2;

/*TODO:

*/


intrinsic MinimalIntermediateIdeals(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]
{Given fractional S-ideals J subset I, returns the minimal (with respect to inclusion) fractional S-ideals K such that J subset K subset I.}
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

intrinsic IntermediateIdeals(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]
{Given fractional S-ideals J subset I, returns all the fractional S-ideals K such that J subset K subset I. They are produced recursively using from the minimal ones.}
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

intrinsic IntermediateIdealsWithPrescribedMultiplicatorRing(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]
{Given fractional S-ideals J subset I, returns all the fractional S-ideals K such that (K:K)=S and  J subset K subset I. They are produced recursively using from the minimal ones.}
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

intrinsic MaximalIntermediateIdeals(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]
{Given fractional S-ideals J subset I, returns the maximal (with respect to inclusion) fractional S-ideals K such that J subset K subset I.}
    assert2 J subset I; // "the ideal J needs to be inside I";
    S:=Order(I);
    assert2 S eq Order(J); // "The ideals must be over the same order";
    if J eq I then 
        return {@ @};
    else
        // If I c K c I is maximal then I/K is a simple module, that is, isomorphic to S/P for some prime P of S.
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

// // Using TraceDualIdeal. Slightly slower 
// intrinsic MaximalIntermediateIdeals(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]
// { Given fractional S-ideals J subset I, returns the maximal (with respect to inclusion) fractional S-ideals K such that J subset K subset I. }
//     assert2 J subset I; // "the ideal J needs to be inside I";
//     S:=Order(I);
//     assert2 S eq Order(J); // "The ideals must be over the same order";
//     if J eq I then 
//         return {@ @};
//     else
//         // I c K c J is a maximal intermediate module iff K^t is a minimal intermediate module J^t c K^t c I^t.
//         min_duals:=MinimalIntermediateIdeals(TraceDualIdeal(J), TraceDualIdeal(I));
//         max:={@ TraceDualIdeal(K) : K in min_duals @};
//         return max;
//     end if;
// end intrinsic;

intrinsic IntermediateIdealsWithTrivialExtension(I::AlgEtQIdl,J::AlgEtQIdl, O::AlgEtQOrd)->SetIndx[AlgEtQIdl]
{Given fractional S-ideals J subset I, returns all the fractional S-ideals K J subset K subset I, and O!!K = I. Note that we need O subset (J:J). They are produced recursively using from the maximal ones.}
    require J subset I : "The ideal J needs to be inside I";
    S:=Order(I);
    require S eq Order(J) : "The ideals must be over the same order";
    require S subset O : "O is not an overorder of Order(I)";
    IO:=O!!I;
    require S!!IO eq I : "J is not an O-ideal";
    queue:={@ I @};
    output:={@  @};
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

intrinsic IntermediateIdealsWithTrivialExtensionAndPrescribedMultiplicatorRing(I::AlgEtQIdl,J::AlgEtQIdl, O::AlgEtQOrd)->SetIndx[AlgEtQIdl]
{Given fractional S-ideals J subset I, returns all the fractional S-ideals K J subset K subset I, O!!K = I, and (K:K) eq S. Note that we need O subset (J:J). They are produced recursively using from the maximal ones.}
    require J subset I : "The ideal J needs to be inside I";
    S:=Order(I);
    require S eq Order(J) : "The ideals must be over the same order";
    require S subset O : "O is not an overorder of Order(I)";
    IO:=O!!I;
    require S!!IO eq I : "J is not an O-ideal";
    queue:={@ I @};
    output:={@  @};
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

/* // this version of the code seems to be a tiny bit slower

intrinsic MinimalIntermediateIdealsVS(I::AlgEtQIdl,J::AlgEtQIdl : primes:=[])->SetIndx[AlgEtQIdl]
{ Given fractional S-ideals J subset I, returns the minimal with respect to inclusion fractional S-ideals K such that J subset K subset I. }
    assert2 J subset I; // "the ideal J needs to be inside I";
    S:=Order(I);
    assert2 S eq Order(J); // "The ideals must be over the same order";
    if J eq I then 
        return {@ @};
    else
        gens_J_over_S:=Generators(J);
        min_ideals:={@ @};
        // Following code is based on the following Lemma.
        // LEMMA: Let M be a finite S-module (eg I/J). Let N be a minimal submodule of M. 
        // Then there exists a maximal ideal P of S such that PN=0;
        // In particular, N is a submodule of the S/P vector space M[P]={ m in M | Pm=0 }.
        // PROOF: PN is a submodule of N. Since N is minimal PN=N or PN=0. 
        // The second case occurs iff P contains Ann_S(N).
        // QED
        ccJImeetS:=ColonIdeal(J,I) meet OneIdeal(S);
        if #primes eq 0 then
            pp:=PrimesAbove(ccJImeetS); 
        else
            pp:=[ P : P in primes | ccJImeetS subset P ];
        end if;
        for P in pp do
            V,v:=QuotientVS(I meet ColonIdeal(J,P),J,P);
            min_Rmod:=MinimalSubmodules(V); //1-dimensional vector spaces
            min_ideals join:={@ Ideal(S, gens_J_over_S cat [ (V!w)@@v : w in Basis(W)]) : W in min_Rmod @};
        end for;
        return min_ideals;
    end if;
end intrinsic;

intrinsic IntermediateIdealsVS(I::AlgEtQIdl,J::AlgEtQIdl)->SetIndx[AlgEtQIdl]
{ Given fractional S-ideals J subset I, returns all the fractional S-ideals K such that J subset K subset I. They are produced recursively using from the minimal ones }
    S:=Order(I);
    require J subset I : "The ideal J needs to be inside I";
    // require S eq Order(J) : "The ideals must be over the same order"; //done in ColonIdeal
    primes:=PrimesAbove(ColonIdeal(J,I) meet OneIdeal(S)); 
    queue:={@ J @};
    output:={@ J @};
    done:={@ @};
    while #queue gt 0 do
        //output join:=queue; //I think this is redundant
        pot_new:=&join[MinimalIntermediateIdealsVS(I,elt : primes:=primes ) : elt in queue ];
        output join:=pot_new;
        done join:=queue;
        queue := pot_new diff done;
    end while;
    return output;
end intrinsic;
*/


/* TESTS

    printf "### Testing IntermediateModules:";
	SetAssertions(2);
	_<x>:=PolynomialRing(Integers());
    f:=x^3-100*x^2-100*x-100;
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    O:=MaximalOrder(K);
    ff:=Conductor(E);
    time _:=MinimalIntermediateIdeals(E!!OneIdeal(O),ff);
    time _:=IntermediateIdeals(E!!OneIdeal(O),ff);
    time _:=IntermediateIdealsWithPrescribedMultiplicatorRing(E!!OneIdeal(O),ff);
    time _:=MaximalIntermediateIdeals(E!!OneIdeal(O),ff);
    time _:=IntermediateIdealsWithTrivialExtension(E!!OneIdeal(O),ff,O);
    time _:=IntermediateIdealsWithTrivialExtensionAndPrescribedMultiplicatorRing(E!!OneIdeal(O),ff,O);
    printf ".";
    f:=x^4+291*x^3-988*x^2-1000*x-1000;
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    O:=MaximalOrder(K);
    ff:=Conductor(E);
    time _:=MinimalIntermediateIdeals(E!!OneIdeal(O),ff);
    time _:=MaximalIntermediateIdeals(E!!OneIdeal(O),ff);
    time _:=IntermediateIdealsWithTrivialExtension(E!!OneIdeal(O),ff,O);
    time _:=IntermediateIdealsWithTrivialExtensionAndPrescribedMultiplicatorRing(E!!OneIdeal(O),ff,O);
    printf ".";
    SetAssertions(1);
    printf " all good!\n"; 

*/
