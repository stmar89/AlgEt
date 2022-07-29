/* vim: set syntax=magma :*/

freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose IntermediateIdeals, 2;

/*TODO:

*/


intrinsic MinimalIntermediateIdeals(I::AlgEtIdl,J::AlgEtIdl)->SetIndx[AlgEtIdl]
{ Given fractional S-ideals J subset I, returns the minimal with respect to inclusion fractional S-ideals K such that J subset K subset I. }
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

intrinsic IntermediateIdeals(I::AlgEtIdl,J::AlgEtIdl)->SetIndx[AlgEtIdl]
{ Given fractional S-ideals J subset I, returns all the fractional S-ideals K such that J subset K subset I. They are produced recursively using from the minimal ones }
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
/* // this version of the code seems to be a tiny bit slower

intrinsic MinimalIntermediateIdealsVS(I::AlgEtIdl,J::AlgEtIdl : primes:=[])->SetIndx[AlgEtIdl]
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

intrinsic IntermediateIdealsVS(I::AlgEtIdl,J::AlgEtIdl)->SetIndx[AlgEtIdl]
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


/* TEST

   

*/
