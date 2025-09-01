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
// Copyright 2025, Stefano Marseglia
/////////////////////////////////////////////////////


freeze;

declare verbose FactPrimes, 3;

declare attributes AlgEtQOrd : SingularPrimes;

///# Primes, factorization and completions

///## Primes
/// Let $S$ be an order in an étale algebra $A$. A `prime` of $S$ is a maximal ideal of $S$, that is, an integral fractional ideal $\mathfrak{p}$ such that $S/\mathfrak{p}$ is a field. 
/// It is well-known that a prime $\mathfrak{p}$ of $S$ is invertible if and only if $\mathfrak{p}$ does not contain the conductor of $S$. Such primes are also called `regular`. A prime which in non-invertible is also called `singular`.
/// In particular, there is only finitely many singular primes.

intrinsic IsPrime(I::AlgEtQIdl) -> BoolElt
{Returns whether the given integral fractional ideal is a prime of its order of definition.}
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

///ditto
intrinsic IsMaximal(I::AlgEtQIdl) -> BoolElt
{Returns whether the given integral fractional ideal is a prime of its order of definition.}
    return IsPrime(I);
end intrinsic;

///ditto
intrinsic IsMaximalIdeal(I::AlgEtQIdl) -> BoolElt
{Returns whether the given integral fractional ideal is a prime of its order of definition.}
    return IsPrime(I);
end intrinsic;

/// Given an integral fractional $S$-ideal $I$, returns the sequence of primes of $S$ containing $I$.
intrinsic PrimesAbove(I::AlgEtQIdl) -> SeqEnum[AlgEtQOrdIdl]
{Given an integral fractional S-ideal I, returns the sequence of primes of S containing I.}
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

intrinsic PlacesAboveRationalPrime(E::AlgEtQ,p::RngIntElt)->SeqEnum[AlgEtQIdl]
{Given an étale algebra and a rational prime, returns the primes of maximal order of the algebra containing the rational prime.}
    if not assigned E`PlacesAboveRationalPrime then
        E`PlacesAboveRationalPrime:=AssociativeArray();
    end if;
    if not IsDefined(E`PlacesAboveRationalPrime,p) then
        require IsPrime(p) : "The integer p needs to be a prime number";
        E`PlacesAboveRationalPrime:=PrimesAbove(p*MaximalOrder(E));
    end if;
    return E`PlacesAboveRationalPrime;
end intrinsic;

intrinsic SingularPrimes(R::AlgEtQOrd) -> SeqEnum[AlgEtQOrdIdl]
{Returns the non-invertible primes of the order.}
    if not assigned R`SingularPrimes then
        R`SingularPrimes:=PrimesAbove(Conductor(R));
    end if;
    return R`SingularPrimes;
end intrinsic;

///ditto
intrinsic NonInvertiblePrimes(R::AlgEtQOrd) -> SetIndx
{Returns the non-invertible primes of the order.}
    return SingularPrimes(R);
end intrinsic;

///## Factorization
/// Every integral fractional $S$-ideal which is coprime to the conductor of $S$ can be written as a product of invertible primes of $S$.
/// This factorization is unique up to permutation of the factors.
/// 
/// Consider the decomposition $K_1\times \cdots \times K_n$ of $A$ into the direct product of its components.
/// Then every prime $\mathfrak{P}$ of the maximal order $\mathcal{O}_A$ of $A$ is generated by a maximal ideal in exactly one component and by $1$ in every other component.
/// It follows that the factorization problem for fractional ideals over the maximal order can be immediately reduced to the number field case.
/// We deal the case of non-maximal orders by taking intersections.

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

/// Given an integral fractional $S$-ideal $I$ coprime with the conductor of $S$ (hence invertible in $S$), returns its factorization into a product of primes of $S$.
intrinsic Factorization(I::AlgEtQIdl) -> Tup
{Given an integral fractional S-ideal I coprime with the conductor of S (hence invertible in S), returns its factorization into a product of primes of S.}
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
    printf " all good!"; 

*/