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

declare attributes AlgEtQOrd:NonGorensteinPrimes, CohenMacaulayType;

///# Gorenstein, Bass and Cohen-Macaulay
/// Let $R$ be an order in an étale algebra $A$ over $\mathbb{Q}$.
/// The `(Cohen-Macaulay) type` of $R$ `at a prime` $\mathfrak{p}$ is defined as minimal number of generators at $\mathfrak{p}$ of the trace dual ideal $R^t$, that is, the dimension $\dim_{R/\mathfrak{p}}(R^t/\mathfrak{p}R^t$.
/// The `(Cohen-Macaulay) type` of $R$ is the maximum of all the local types.
///  
/// If the type of $R$ at $\mathfrak{p}$ is one, we say that $R$ is `Gorenstein at the prime` $\mathfrak{p}$.
/// It follows that if $\mathfrak{p}$ is invertible then $R$ is Gorenstein at $\mathfrak{p}$.
/// Moreover, $R$ is Gorenstein at $\mathfrak{p}$ if and only if every ideal with multiplicator ring $R$ is locally principal at $\mathfrak{p}$.
/// The order $R$ is said to be `Gorenstein` if it is so at all its (singular) primes.
/// This is equivalent to every fractional $R$-ideal $I$ with $(I:I)=R$ being invertible.
///  
/// The order $R$ is `Bass at the prime` $\mathfrak{p}$ if every overorder $S$ of $R$ is Gorenstein at the finitely many primes above $\mathfrak{p}$. The order $R$ is `Bass` if it is so at all (singular) primes.
/// One can show that $R$ is Bass at $\mathfrak{p}$ if and only if $\dim_{R/\mathfrak{p}}(\mathcal{O_A}/\mathfrak{p}\mathcal{O}_A)\leq 2$, where $\mathcal{O}_A$ is the maximal order of $A$.
/// 
/// Reference: Stefano Marseglia, "Cohen-Macaulay type of orders, generators and ideal classes" (Journal of Algebra 658 (2024), 247-276.)

//------------
// CohenMacaulayType
//------------

/// Given an order $S$ and a prime $P$, returns the Cohen-Macaulay type of $S$ at $P$.
intrinsic CohenMacaulayTypeAtPrime(S::AlgEtQOrd,P::AlgEtQIdl)->RngIntElt
{Given an order S and a prime P, returns the Cohen-Macaulay type of S at P. This integer equals the dimension of S^t/P*S^t where S^t is the trace dual of S.}
    if IsGorenstein(S) then
        return 1;
    end if;
    if assigned S`NonGorensteinPrimes then
        pp,dPs:=NonGorensteinPrimes(S);
        i:=Index(pp,P);
        if i eq 0 then // S is Gorenstein at P
            return 1;
        else
            return dPs[i];
        end if;
    end if;
    // No early exit ocured. We compute it.
    St:=TraceDualIdeal(S);
    k:=Integers() ! Index(S,P);
    v:=Integers() ! Index(St,P*St);
    dP:=Ilog(k,v);
    return dP;
end intrinsic;

/// Given an order $S$, returns its Cohen-Macaulay type.
intrinsic CohenMacaulayType(S::AlgEtQOrd)->RngIntElt
{Given an order S, returns its Cohen-Macaulay type. This integer equals the max dimension of S^t/P*S^t where S^t is the trace dual of S and P runs over all (non-Gorenstein) primes of S.}
    if not assigned S`CohenMacaulayType then
        pp,dps:=NonGorensteinPrimes(S);
        if #pp eq 0 then
            S`CohenMacaulayType:=1; //S is Gorenstein
        else
            S`CohenMacaulayType:=Max(dps);
        end if;
    end if;
    return S`CohenMacaulayType;
end intrinsic;

//------------
// Gorenstein
//------------

/// Given an order $S$ and a prime $P$, returns whether $S$ is Gorenstein at $P$.
intrinsic IsGorensteinAtPrime(S::AlgEtQOrd,P::AlgEtQIdl) -> BoolElt
{Given an order S and a prime P, returns whether S is Gorenstein at P, that is, if every fractional ideal I with (I:I)=S is locally principal at P.}
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

/// Returns whether the given order is Gorenstein.
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
end intrinsic;

//------------
// Bass
//------------


/// Given an order $S$ and prime $P$ of $S$, returns whether $S$ is Bass at $P$.
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

/// Returns whether the given order is Bass.
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


//------------
// NonGorensteinPrimes
//------------

/// Given an order $S$, returns two sequences: the first containis the primes at which $S$ is locally not Gorenstein; the second containis the Cohen-Macaulay types of $S$ at these primes.
intrinsic NonGorensteinPrimes(S::AlgEtQOrd)->SeqEnum,SeqEnum
{Given an order S, returns two sequences: the first containis the primes at which S is locally not Gorenstein; the second containis the CohenMacaulay types of S at these primes, that is, the dimension of S^t/PS^t over S/P, where S^t is the TraceDualIdeal of S.}
    if not assigned S`NonGorensteinPrimes then
        if IsGorenstein(S) then
            S`NonGorensteinPrimes:=<[],[]>;
        else
            St:=TraceDualIdeal(S);
            pp0:=PrimesAbove(St*ColonIdeal(S,St));
            pp:=[];
            dPs:=[];
            for P in pp0 do
                k:=Integers() ! Index(S,P);
                v:=Integers() ! Index(St,P*St);
                dP:=Ilog(k,v);
                assert dP gt 1; // don't want P such that S_P is Gorenstein
                Append(~pp,P);
                Append(~dPs,dP);
            end for;
            assert Min(dPs) ne 1;
            S`NonGorensteinPrimes:=< pp,dPs >;
        end if;
    end if;
    return Explode(S`NonGorensteinPrimes);
end intrinsic;

///# Example 6
/// ```
/// // Given an étale algebra A of dimension n over Q, the maximum Cohen-Macaulay type that an order can have is n-1.
/// // An order with such maximal Cohen-Macaualay type can be found among the overorders of the order Z+p*O_A where p an arbitrary rational prime. We verify this statement in an example.
/// _<x>:=PolynomialRing(Integers());
/// f:=x^4+8;
/// A:=EtaleAlgebra(f);
/// O:=MaximalOrder(A);
/// p:=7;
/// E:=Order([p*z:z in ZBasis(O)]);
/// n:=Dimension(A);
/// oo:=OverOrders(E);
/// exists{S:S in oo|CohenMacaulayType(S) eq n-1};
/// ```


/////////////////////////////////////////////////////
// functions for orders that locally have CM-type \leq 2
/////////////////////////////////////////////////////
// if S^t/P*S^t has dimension 2 over S/P then,
// locally at P, the only fractional ideals with 
// multiplicator ring S are S and S^t.

wkicm_bar_CM_type2:=function(S,pp)
// Input:
//  S is an order such that for all primes P we have d(P) le 2.
//  St is the trace dual of S
//  pp is the list of primes P at which d(P)=2.
// Output: a sequence of ideals of S containing all weak equivalence 
//  classes of fractional ideals with multiplicator ring S.
    St:=TraceDualIdeal(S);
    if #pp eq 1 then
        return [OneIdeal(S),St];
    else
        I:=MakeIntegral(St);
        pows:=[];
        for P in pp do
            _,p,e:=IsPrimePower(Integers()!Index(S,P)); // e=Valuation(p,Index(S,P));
            m:=Valuation(Index(S,I),p) div e; // Here we are using Proposition 4.2 of Klueners and Pauli 
                                              // "Computing residue class rings ...".
                                              // With such an m we have:
                                              //   pp[i]^m \subseteq I_pp[i]
            Append(~pows,P^m);
        end for;
        // coprime ideals, so we can use meet instead of *, which is faster
        pows_hat:=[ &meet[ pows[j] : j in [1..#pp] | j ne i ] : i in [1..#pp] ];
        Is:=[ seq : seq in CartesianProduct([[OneIdeal(S),I] : i in [1..#pp]]) ];
        out:=[];
        for I in Is do
            L:=&+[ (I[i]+pows[i]) * pows_hat[i] : i in [1..#pows_hat]];
            Append(~out,L);
        end for;
        return out;
     end if;
end function;

/* TESTS

    printf "### Testing LowCohenMacaulayType:";
    SetAssertions(2);
    _<x>:=PolynomialRing(Integers());
    f:=x^4-10000*x^3-10000*x^2-10000*x-10000; 
    //AttachSpec("~/packages_github/AlgEt/spec");
    A:=EtaleAlgebra(f);
    E:=EquationOrder(A);
    oo:=FindOverOrders(E); // ~13 secs
    assert #oo eq 297;
    for S in oo do
        pp:=SingularPrimes(S);
        for P in pp do        
            _:=CohenMacaulayTypeAtPrime(S,P);
        end for;
        pp_ng:=NonGorensteinPrimes(S);
        assert pp_ng subset pp;
        _:=CohenMacaulayType(S);
        printf ".";
    end for;
    SetAssertions(1);
    printf " all good!"; 
*/
