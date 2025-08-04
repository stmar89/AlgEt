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

declare verbose IdealsOfIndex, 2;

intrinsic IdealsOfIndex(O::RngOrd, N::RngIntElt) -> SeqEnum[RngOrdIdl]
{Given an order O in a number field and a positive integer N, returns all the ideals I of index [O:I]=N.}
  vprintf IdealsOfIndex,2 : "IdealsOfIndex RngOrd Int\n";
  require N gt 0 : "N must be a strictly positive integer";
  if N eq 1 then
        return [1*O];
  end if;
  return [ J : J in IdealsUpTo(N, O) | Norm(J) eq N ];
end intrinsic;

intrinsic IdealsOfIndex(I::RngOrdIdl, N::RngIntElt) -> SeqEnum[RngOrdIdl]
{Given an ideal I in an order O in a number field and a positive integer N, with N coprime with the conductor, returns all the ideals J contained in I with index [I:J]=N.}
    vprintf IdealsOfIndex,2 : "IdealsOfIndex RngOrdIdl\n";
    require N gt 0 : "N must be a strictly positive integer";
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
    vprintf IdealsOfIndex,2 : "IdealsOfIndex RngOrdFracIdl\n";
    require N gt 0 : "N must be a strictly positive integer";
    if N eq 1 then
        return [I];
    end if;
    d := Denominator(I);
    dI := Order(I)!!(d*I);
    // coprimlity with the conductor is tested here
    Js := IdealsOfIndex(dI, N);
    return [J/d : J in Js];
end intrinsic;

ideals_of_index_product:=function(Is, N)
    // Given a list Is of ideals representing the ideal I as a product and positive integer N, 
    // returns all the ideals of index N. N has to be coprime with the conductor.
    vprintf IdealsOfIndex,2 : "IdealsOfIndexProduct internal\n";
    if #Is eq 1 then
        return [<elt> : elt in IdealsOfIndex(Is[1], N)];
    end if;
    result := [];
    for d in Divisors(N) do
        J1 := IdealsOfIndex(Is[1], N div d);
        Jk := $$(<Is[i] : i in [2..#Is]>, d);
        if #J1 ge 0 and #Jk ge 0 then
            for c in CartesianProduct(J1, Jk) do
                Append(~result, <c[1]> cat c[2]);
            end for;
        end if;
    end for;
    return result;
end function;

intrinsic IdealsOfIndex(I::AlgEtQIdl, N::RngIntElt : Method := "Default") -> SeqEnum[AlgEtQIdl]
{Given an O-ideal I in O and a positive integer N, returns all the subideals J of I with index [I:J]=N. The function is very fast if N is coprime to the conductor of O. If this conditions are not satisfied a slow algorithm is used which doesn't require additional hypothesis. One can force the slow by setting the vararg Method:="Slow".}
    require Method in {"Default","Slow"} : "The method is not recognized. It should be either Default or Slow"; 
    require N gt 0 : "N must be a strictly positive integer";
    if N eq 1 then
        return [I];
    end if;
    if Method eq "Slow" then
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
            JA:= Order(I)!!Ideal(MultiplicatorRing(I) , J);
            assert Index(I,JA) eq N;
            Append(~result,JA);
        end for;
        return result;
    else
        vprintf IdealsOfIndex,2 : "Slow version\n";
        result:=[ K : K in IntermediateIdealsOfIndex(I,N*I,N)]; // conver to SeqEnum
        return result;
    end if;
end intrinsic;

intrinsic IdealsOfIndex(O::AlgEtQOrd, N::RngIntElt : Method := "Default") -> SeqEnum[AlgEtQIdl]
{Given an order O and a positive integer N, returns all the O-ideals J with index [O:J]=N. The function is very fast if N is coprime to the conductor of O. If this conditions are not satisfied a slow algorithm is used which doesn't require additional hypothesis. One can force the slow-naive by setting the vararg Method:="Slow".}
    return IdealsOfIndex(OneIdeal(O),N : Method:=Method);
end intrinsic;

/* TESTS

    printf "### Testing IdealsOfIndex:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    SetVerbose("IdealsOfIndex",1);
    SetAssertions(2);

    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    E:=EquationOrder(A);
    O:=MaximalOrder(A);
    ind:=Index(O,E);
    for N in [1..15] do
        printf "."; 
        // test with maximal order
        assert Seqset(IdealsOfIndex(O,N)) eq Seqset(IdealsOfIndex(O,N : Method:="Slow"));

        //test with equation order
        if IsCoprime(N,ind) then
            assert Seqset(IdealsOfIndex(E,N)) eq Seqset(IdealsOfIndex(E,N : Method:="Slow"));
        else
            _:=IdealsOfIndex(E,N);
        end if;
    end for;
    SetAssertions(1);    
    printf " all good!"; 
*/
