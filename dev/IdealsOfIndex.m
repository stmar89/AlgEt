/* vim: set syntax=magma :*/

freeze;

/////////////////////////////////////////////////////
// Ideals of given index for Etale Q algebras
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl, http://www.staff.science.uu.nl/~marse004/
// Edgar Costa, MIT
/////////////////////////////////////////////////////

declare verbose IdealsOfIndex, 2;

intrinsic IdealsOfIndex(O::RngOrd, N::RngIntElt) -> SeqEnum[RngOrdIdl]
{Given an Order, retuns all the ideals of index N in that order}
  vprintf IdealsOfIndex : "IdealsOfIndex RngOrd Int\n";
  if N eq 1 then
        return [O];
  end if;
  Js := IdealsUpTo(N, O);
  // Js are ordered by norm, and we only care about the ones with Norm = N * norm_I
  result := [];
  for J in Reverse(Js) do
    if Norm(J) eq N then
        Append(~result, J);
    else
        break;  //the other ideals in Js will have norm < N.
    end if;
  end for;
  return result;
end intrinsic;

intrinsic IdealsOfIndex(I::RngOrdIdl, N::RngIntElt) -> SeqEnum[RngOrdIdl]
{Given an ideal I in an order O in a number field and a positive integer N, with N coprime with the conductor, returns all the ideals J contained in I with index [I:J]=N.}
    vprintf IdealsOfIndex : "IdealsOfIndex RngOrdIdl\n";
    if N eq 1 then
        return [I];
    end if;
    O := Order(I);
    OK := MaximalOrder(NumberField(O));
    index_OK_O := Index(OK, O);
    require N gt 0 and GCD(index_OK_O, N) eq 1 : "N is not coprime with the conductor of Order(I)";
    Js := IdealsOfIndex(OK, N);
    ff:=OK !! Conductor(O);
    assert forall{J : J in Js | J+ff eq 1*OK};
    result := [];
    for J in Js do
            K := (O meet J) * I; // OK/J=O/(J meet O)=I/K, where the second isomorphism holds because 
                                 // (J meet O) is invertible in O, since it is coprime with ff.
            Append(~result, K);
    end for;
    return result;
end intrinsic;

intrinsic IdealsOfIndex(I::RngOrdFracIdl, N::RngIntElt) -> SeqEnum[RngOrdFracIdl]
{Given an ideal I in an order O in a number field and a positive integer N, with N coprime with the conductor, returns all the ideals J contained in I with index [I:J]=N.}
    vprintf IdealsOfIndex : "IdealsOfIndex Frac\n";
    if N eq 1 then
        return [I];
    end if;
    d := Denominator(I);
    dI := Order(I)!!(d*I);
    Js := IdealsOfIndex(dI, N);
    return [J/d : J in Js];
end intrinsic;

ideals_of_index_product:=function(Is, N)
    // Given a list Is of ideals representing the ideal I as a product and positive integer N, 
    // returns all the ideals of index N. N has to be coprime with the conductor.
    vprintf IdealsOfIndex : "IdealsOfIndexProduct\n";
    if #Is eq 1 then
        return [<elt> : elt in IdealsOfIndex(Is[1], N)];
    end if;
    result := [];
    for d in Divisors(N) do
        J1 := IdealsOfIndex(Is[1], Integers()!(N/d));
        Jk := $$(<Is[i] : i in [2..#Is]>, d);
        if #J1 ge 0 and #Jk ge 0 then
            for c in CartesianProduct(J1, Jk) do
                assert #c[2] eq #Jk;
                Append(~result, <c[1]> cat c[2]);
            end for;
        end if;
    end for;
    return result;
end function;

intrinsic IdealsOfIndex(I::AlgEtQIdl, N::RngIntElt : Method := "Default") -> SeqEnum[AlgEtQIdl]
{Given an O-ideal I in O and integer N returns all the subideals J of I with index [I:J]=N. The function is very fast if N is coprime to the conductor of O. If this conditions are not satisfied a slow algorithm is used which doesn't require additional hypothesis. One can force the slow-naive by setting the vararg Method:="Naive".}
    require Method in {"Default","Naive"} : "The method is not recognized. It should be either Default or Naive"; 
    if N eq 1 then
        return [I];
    end if;
    if Method eq "Naive" then
        test := false;
    else
        test, dec := IsProductOfIdeals(MultiplicatorRing(I)!!I);
        for J in dec do
            O := Order(J);
            OK := MaximalOrder(NumberField(O));
            index_OK_O := Index(OK, O);
            if GCD(index_OK_O, N) ne 1 then
                test :=false;
                break;
            end if;
        end for;
    end if;
    if test then
        Js := ideals_of_index_product(dec, N);
        O := Order(I);
        A := Algebra(O);
        result := [];
        for J in Js do
            assert #J eq #dec;
            assert #J eq #Components(A);
            gen_inA := [];
            for i := 1 to #J do
                comps,embs := Components(A);
                gen_inA := gen_inA cat [embs[i](y) : y in Basis(J[i], comps[i])];
            end for;
            JA:= Ideal(Order(I) , gen_inA);
            assert Index(I,JA) eq N;
            Append(~result,JA);
        end for;
        return result;
    else
        vprintf IdealsOfIndex : "Naive version\n";
        result:=[ K : K in IntermediateIdealsOfIndex(I,N*I,N)]; // converto SeqEnum
        return result;
    end if;
end intrinsic;

// MOVE ME TO IntermediateIdeals.m
intrinsic IntermediateIdealsOfIndex(I::AlgEtQIdl,J::AlgEtQIdl,N::RngIntElt)->SetIndx[AlgEtQIdl]
{Given ideals J subset I over the same order, and a positive integer N, it returns all the ideals K such that J subset K subset I and [I:K]=N. Theser are produced by recursively searching for maximal submodules.}
    require J subset I : "The ideal J needs to be inside I";
    require N gt 0 : "N must be a strictly positive integer";
    S:=Order(I);
    require S eq Order(J) : "The ideals must be over the same order";
    // early exits
    if Index(I,J) eq N then
        return {@ J @};
    elif (Index(I,J) mod N) ne 0 then
        output:={@ Universe({@ I @}) | @}; //empty set
        return output;
    end if;
    // we start the recursion
    queue:={@ I @};
    output:={@  @};
    done:={@ @};
    while #queue gt 0 do
        pot_new:=&join[MaximalIntermediateIdeals(elt,J) : elt in queue ];
        output join:={@ K : K in pot_new | Index(I,K) eq N @};
        done join:=queue;
        queue := {@ M : M in pot_new diff done | Index(I,M) lt N @}; // if [I:M] ge N then for all K c M 
                                                                   // we have that [I:K]>N. We don't want such M's 
                                                                   // in the queue.
    end while;
    return output;
end intrinsic;

