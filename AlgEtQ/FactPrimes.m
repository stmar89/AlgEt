/* vim: set syntax=magma :*/

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

declare verbose FactPrimes, 3;

declare attributes AlgEtQOrd : SingularPrimes;

//----------
// Factorization and Prime
//----------

function factorizationMaximalOrder(I)
//given an ideal of the maximal order of an étale algebra, returns the factorization into a product of prime ideals
    O:=Order(I);
    assert IsMaximal(O);
    test,OasProd:=IsProductOfOrders(O);
    assert test; //since we assume that O is maximal
    A:=Algebra(O);
    test,IasProd:=IsProductOfIdeals(I);
    assert test;
    fac:=[]; //this will be the factorization of I
    nf,embs:=Components(A);
    tup_one_ideals:=< 1*O : O in OasProd >;
    for i in [1..#nf] do
        IL:=IasProd[i];
        if not One(Order(IL)) in IL then
            facL:=Factorization(IL); // < <P,e> : ... >;
            assert #facL gt 0;
            for p in facL do
                tup_p:=tup_one_ideals;
                tup_p[i]:=p[1]; // we replace the i-th ideal with P
                assert2 tup_p ne tup_one_ideals;
                P:=Ideal(O,tup_p);
                ZBasisLLL(P);
                assert2 P ne OneIdeal(O);
                P`IsPrime:=true; //we know P is prime
                assert2 IsPrimePower(Integers() ! Index(O,P));
                Append(~fac,<P,p[2]>);
            end for;
         end if;
    end for;
    assert2 I eq &*[p[1]^p[2] : p in fac];
    return fac;
end function;

intrinsic Factorization(I::AlgEtQIdl) -> Tup
{Given an integral S-ideal I coprime with the conductor of S (hence invertible in S), returns its factorization into a product of primes of S.}
    S:=Order(I);
    require IsIntegral(I) and I ne OneIdeal(S): "the argument must be a proper integral ideal";
    if not assigned I`Factorization then    
        if IsMaximal(S) then
            I`Factorization:=factorizationMaximalOrder(I);
        else
            fS:=Conductor(S);
            require IsCoprime(fS,I): "the ideal must be coprime with the conductor of the order of definition";
            O:=MaximalOrder(Algebra(I));
            IO:=O !! I;
            facO:=factorizationMaximalOrder(IO);
            primesO:=[ p[1] : p in facO ];
            primesS:=Setseq({ OneIdeal(S) meet (S!!PO) : PO in primesO }); //cancel the doubles
            facS:=<>;
            for i in [1..#primesS] do
                P:=primesS[i];
                ZBasisLLL(P);
                P`IsPrime:=true;
                expP:=&+([ pO[2] : pO in facO | (S meet (S!!pO[1])) eq P ]);
                Append(~facS, <P,expP>);
            end for;
            assert2 I eq &*([ p[1]^p[2] : p in facS ]);
            I`Factorization:=facS;
        end if;
     end if;
     return I`Factorization;
end intrinsic;

intrinsic PrimesAbove(I::AlgEtQIdl) -> SeqEnum[AlgAssEtOrdIdl]
{Given an integral S-ideal, returns the sequence of maximal ideals P of S above I.}
    require IsIntegral(I): "the ideal must be integral";
    if not assigned I`PrimesAbove then
        S:=Order(I);
        if assigned I`Factorization then
            primes:=[P[1] : P in I`Factorization];
        elif I eq OneIdeal(S) then
            primes:=[];
        else
            if IsMaximal(S) then
                O:=S;
                IO:=I;
                fac:=Factorization(IO);
                primes:=[P[1] : P in fac]; //they are all distinct
            else
                O:=MaximalOrder(Algebra(I));
                IO:=O!!I;
                fac:=Factorization(IO);
                primes:= Setseq({ OneIdeal(S) meet (S!!PO[1]) : PO in fac }); //remove doubles
                for i in [1..#primes] do
                    P:=primes[i];
                    P`IsPrime:=true;
                end for;
            end if;
            assert2 forall{P : P in primes | IsIntegral(P)};
            assert2 forall{P : P in primes | I subset P};
            assert2 forall{P : P in primes | IsPrimePower(Integers() ! Index(S,P))};
        end if;
        for P in primes do
            ZBasisLLL(P);
        end for;
        I`PrimesAbove:=primes;
    end if;
    return I`PrimesAbove;
end intrinsic;

intrinsic SingularPrimes(R::AlgEtQOrd) -> SeqEnum[AlgAssEtOrdIdl]
{Returns the non-invertible primes of the order.}
    if not assigned R`SingularPrimes then
        R`SingularPrimes:=PrimesAbove(Conductor(R));
    end if;
    return R`SingularPrimes;
end intrinsic;

intrinsic NonInvertiblePrimes(R::AlgEtQOrd) -> SetIndx
{Returns the non-invertible primes of the order.}
    return SingularPrimes(R);
end intrinsic;

intrinsic IsPrime(I::AlgEtQIdl) -> BoolElt
{Given an integral S-ideal, returns if the ideal is a prime fractional ideal of S, that is a maximal S ideal.}
    require IsIntegral(I): "the ideal must be integral";
    if not assigned I`IsPrime then
        prim:=PrimesAbove(I);
        if #prim eq 1 and I eq prim[1] then
            bool:=true;
        else
            bool:=false;
        end if;
        I`IsPrime:=bool;
    end if;
    return I`IsPrime;
end intrinsic;

intrinsic IsBassAtPrime(S::AlgEtQOrd,P::AlgEtQIdl) -> BoolElt
{Check if the order is Bass at the prime ideal P, that is, if every overorder of S is Gorenstein at the primes above P.}
// we compute the maximal order O and check if O/PO is at most 2-dimensional over S/P
// see Hofman,Sircana, On the computations of over-orders, Definition 5.23
    if assigned S`IsBass and S`IsBass then
        return true;
    end if;
    if IsMaximal(S) then // we need to compute the maximal order in any case.
        return true;
    else
        O:=MaximalOrder(Algebra(S));
        k:=Integers() ! Index(S,P);
        OS:=S!!OneIdeal(O);
        N:=Integers() ! Index(OS,P*OS);
        //N = k^(dim_P)
        assert N mod k eq 0;
        dim_P:=Ilog(k,N);
        return dim_P le 2;
    end if;
end intrinsic;

intrinsic IsBass(S::AlgEtQOrd) -> BoolElt
{Check if the order is Bass, that is, if every overorder of S is Gorenstein.}
// we compute the maximal order O and check if O/PO is at most 2-dimensional over S/P for every singular prime P
// This coincides with the usual definition since O has the maximal number of generators as a fractional S ideal and S is Bass iff every ideal can be generated by at most 2-elements.
// see Hofman,Sircana, On the computations of over-orders, Definition 5.23
    if not assigned S`IsBass then
        if IsMaximal(S) then 
            S`IsBass:=true;
        else
            O:=MaximalOrder(Algebra(S));
            sing:=SingularPrimes(S);
            S`IsBass:=forall{P:P in sing | IsBassAtPrime(S,P)};
        end if;
    end if;
    return S`IsBass;
end intrinsic;

intrinsic IsGorensteinAtPrime(S::AlgEtQOrd,P::AlgEtQIdl) -> BoolElt
{Check if the order is Gorenstein at the prime ideal P, that is, if every fractional ideal I with (I:I)=S is locally principal at P.}
// Let St be the trace dual ideal of S. We check if St/PSt is one dimensional over S/P
    if assigned S`IsGorenstein and S`IsGorenstein then
        return true;
    end if;
    if assigned S`IsMaximal and S`IsMaximal then //we dont want to trigger the computation of the maximal order
        return true;
    else
        k:=Integers() ! Index(S,P);
        St:=TraceDualIdeal(S);
        N:=Integers() ! Index(St,P*St);
        //N = k^(dim_P)
        assert N mod k eq 0;
        dim_P:=Ilog(k,N);
        return dim_P eq 1;
    end if;
end intrinsic;


intrinsic IsGorenstein(O::AlgEtQOrd)->BoolElt
{Checks if the order O is Gorenstein, that is if the TraceDualIdeal of O is invertible, or equivalently, if all fractional ideals I with (I:I)=O are invertible.}
    if not assigned O`IsGorenstein then
        if assigned O`IsMaximal and O`IsMaximal then
            O`IsGorenstein:=true;
        else
            T:=TraceDualIdeal(O);
            O`IsGorenstein:=IsInvertible(T);
        end if;
    end if;
    return O`IsGorenstein;
end intrinsic

/* TESTS

    printf "### Testing Primes and Factorizaton:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    SetAssertions(2);

    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16);
    A:=EtaleAlgebra(f);
    E1:=EquationOrder(A);
    ff:=Conductor(E1);
    assert PrimesAbove(Conductor(E1)) eq SingularPrimes(E1);
    printf ".";

    f:=x^6 + 8*x^5 + 50*x^4 + 200*x^3 + 1250*x^2 + 5000*x + 15625;
    A:=EtaleAlgebra(f);
    assert #SingularPrimes(MaximalOrder(A)) eq 0;
    R:=EquationOrder(A);
    assert #SingularPrimes(R) eq 5;
    assert #NonInvertiblePrimes(R) eq 5;
    printf ".";

    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    E1:=EquationOrder(A);
    E2:=ProductOfEquationOrders(A);
    
    _:=PrimesAbove(Conductor(E1));
    _:=PrimesAbove(Conductor(E2));
    assert IsGorenstein(E1);
    assert IsGorenstein(E2);

    ids:=[ Ideal(E1,[Random(E1) : i in [1..10]]) : i in [1..100]];
    ids0:=[ I : I in ids | I ne OneIdeal(E1) and IsInvertible(I) ];
    ids:=[];
    for I in ids0 do
        _,J:=CoprimeRepresentative(I,Conductor(E1));
        Append(~ids,J);
    end for;
    facs:=[ Factorization(I) : I in ids ];
    printf ".";
    SetAssertions(1);
    printf " all good!\n"; 

*/
