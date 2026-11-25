/////////////////////////////////////////////////////
// Copyright 2025.
// Stefano Marseglia, stefano.marseglia89@gmail.com
// https://stmar89.github.io/index.html
//
// Distributed under the terms of the CC-BY 4.0 licence.
// https://creativecommons.org/licenses/by/4.0/
/////////////////////////////////////////////////////



freeze;

declare verbose FactPrimes, 3;

declare attributes AlgEtQOrd : SingularPrimes,
                               PrimesBelow; // The cache of PrimeBelow: an AssociativeArray keyed by
                                            // the orders T
                                            // containing S, whose value at T is an AssociativeArray
                                            // keyed by the primes P of T, with value
                                            // PrimeBelow(S,P). The primes are grouped by their
                                            // order because ideals over different orders cannot be
                                            // compared.

declare attributes AlgEtQ : PlacesAboveRationalPrime;

///# Primes, factorization and completions

///## Primes
/// Let $S$ be an order in an étale algebra $A$. A `prime` of $S$ is a maximal ideal of $S$, that is, an integral fractional ideal $\mathfrak{p}$ such that $S/\mathfrak{p}$ is a field.
/// It is well-known that a prime $\mathfrak{p}$ of $S$ is invertible if and only if $\mathfrak{p}$ does not contain the conductor of $S$. Such primes are also called `regular`. A prime which is non-invertible is also called `singular`.
/// In particular, there are only finitely many singular primes.

intrinsic IsPrime(I::AlgEtQIdl) -> BoolElt
{Returns whether the given integral fractional ideal is a prime of its order of definition.}
    if not assigned I`IsPrime then
        require IsIntegral(I): "the ideal must be integral";
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
{ " } // "
    return IsPrime(I);
end intrinsic;

///ditto
intrinsic IsMaximalIdeal(I::AlgEtQIdl) -> BoolElt
{ " } // "
    return IsPrime(I);
end intrinsic;

/// Given an order $S$ and a prime $P$ of an order $T$ containing $S$, returns the prime $S \cap P$ of $S$ below $P$. The output is stored in $S$: asking again for the same prime returns the same ideal, with all its attributes.
intrinsic PrimeBelow(S::AlgEtQOrd, P::AlgEtQIdl) -> AlgEtQIdl
{Given an order S and a prime P of an order T containing S, returns the prime S meet P of S below P. The output is stored in S: asking again for the same prime returns the same ideal, with all its attributes.}
    T:=Order(P);
    if T eq S then
        return P;
    end if;
    if assigned S`PrimesBelow and IsDefined(S`PrimesBelow,T) and IsDefined(S`PrimesBelow[T],P) then
        return S`PrimesBelow[T][P];
    end if;
    require IsPrime(P) : "the ideal must be prime";
    // The inclusion test is skipped when T is the maximal order, which contains every order: this
    // is the main use, and with many orders S the test would cost more than the intersection
    // itself.
    require (assigned T`IsMaximal and T`IsMaximal) or S subset T : "the order must be contained in the order of the prime";
    // S/(S meet P) embeds in the finite field T/P, so it is a finite domain, hence a field: S meet
    // P is a prime of S.
    Q:=OneIdeal(S) meet S!!P;
    Q`IsPrime:=true;
    cache:=assigned S`PrimesBelow select S`PrimesBelow else AssociativeArray();
    // Distinct primes, also of different orders, can be above the same prime of S: if it is already
    // stored, we return it, so that each prime of S is one ideal, whose attributes (residue field,
    // ...) are computed once.
    for T2 in Keys(cache) do
        for P2 in Keys(cache[T2]) do
            if cache[T2][P2] eq Q then
                Q:=cache[T2][P2];
                break T2;
            end if;
        end for;
    end for;
    below:=IsDefined(cache,T) select cache[T] else AssociativeArray();
    below[P]:=Q;
    cache[T]:=below;
    S`PrimesBelow:=cache;
    return Q;
end intrinsic;

/// Given an integral fractional $S$-ideal $I$, returns the sequence of primes of $S$ containing $I$.
intrinsic PrimesAbove(I::AlgEtQIdl) -> SeqEnum[AlgEtQIdl]
{Given an integral fractional S-ideal I, returns the sequence of primes of S containing I.}
    if not assigned I`PrimesAbove then
        require IsIntegral(I): "the ideal must be integral";
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
                // The set removes the repetitions (several primes of O can be above the same prime
                // of S). Its order depends on the hash values of the primes, so the order of the
                // output of PrimesAbove is not otherwise meaningful.
                primes:= Setseq({ PrimeBelow(S,PO[1]) : PO in fac });
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
        E`PlacesAboveRationalPrime[p]:=PrimesAbove(p*MaximalOrder(E));
    end if;
    return E`PlacesAboveRationalPrime[p];
end intrinsic;

intrinsic SingularPrimes(R::AlgEtQOrd) -> SeqEnum[AlgEtQIdl]
{Returns the non-invertible primes of the order.}
    if not assigned R`SingularPrimes then
        R`SingularPrimes:=PrimesAbove(Conductor(R));
    end if;
    return R`SingularPrimes;
end intrinsic;

///ditto
intrinsic NonInvertiblePrimes(R::AlgEtQOrd) -> SeqEnum[AlgEtQIdl]
{ " } // "
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
//given an ideal of the maximal order of an étale algebra, returns the factorization into a product
//of prime ideals
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
    // The powers are formed by repeated products, not with '^': for an invertible ideal '^' calls
    // TwoGeneratingSet, which can call CoprimeRepresentative, PrimesAbove and Factorization of an
    // equal ideal again, an infinite recursion at assertion level 2 (Magma then crashes).
    assert2 I eq &*[ &*[ p[1] : i in [1..p[2]] ] : p in fac ];
    return fac;
end function;

/// Given an integral fractional $S$-ideal $I$ coprime with the conductor of $S$ (hence invertible in $S$), returns its factorization into a product of primes of $S$.
intrinsic Factorization(I::AlgEtQIdl) -> Tup
{Given an integral fractional S-ideal I coprime with the conductor of S (hence invertible in S), returns its factorization into a product of primes of S.}
    if not assigned I`Factorization then
        S:=Order(I);
        require IsIntegral(I) and I ne OneIdeal(S): "the argument must be a proper integral ideal";
        if IsMaximal(S) then
            I`Factorization:=factorizationMaximalOrder(I);
        else
            fS:=Conductor(S);
            require IsCoprime(fS,I): "the ideal must be coprime with the conductor of the order of definition";
            O:=MaximalOrder(Algebra(I));
            IO:=O !! I;
            facO:=factorizationMaximalOrder(IO);
            // I is coprime with the conductor of S. So every prime P of S containing I is below
            // exactly one prime PO of O containing IO, and the valuations of I at P and of IO at PO
            // are equal, since S and O are equal locally at P. Hence we only need to intersect each
            // PO with S.
            facS:=<>;
            for pO in facO do
                P:=PrimeBelow(S,pO[1]);
                ZBasisLLL(P);
                P`IsPrime:=true;
                Append(~facS, <P,pO[2]>);
            end for;
            assert2 #{ p[1] : p in facS } eq #facS; // the primes of S are distinct
            assert2 I eq &*[ &*[ p[1] : i in [1..p[2]] ] : p in facS ]; // no '^': see factorizationMaximalOrder
            I`Factorization:=facS;
        end if;
     end if;
     return I`Factorization;
end intrinsic;


/* TESTS

    printf "### Testing Primes and Factorizaton:";
    //AttachSpec("../spec");
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

    // random principal ideals: an ideal generated by many random elements is almost always E1
    ids:=[ Random(E1)*E1 : i in [1..5] ];
    ids0:=[ I : I in ids | I ne OneIdeal(E1) and IsInvertible(I) ];
    ids:=[];
    for I in ids0 do
        _,J:=CoprimeRepresentative(I,Conductor(E1));
        // J can be E1 itself, for example when the inverse of a generator of I is short
        if J ne OneIdeal(E1) then
            Append(~ids,J);
        end if;
    end for;
    facs:=[* Factorization(I) : I in ids *]; // a list: tuples of different lengths are not in the
                                             // same universe
    // an ideal with many prime factors, coprime with the conductor.
    // We compare its primes with PrimesAbove of a copy without the attribute Factorization.
    ps:=[ p : p in PrimesUpTo(20) | IsCoprime(p*OneIdeal(E1),Conductor(E1)) ];
    I:=(&*ps)*OneIdeal(E1);
    fac:=Factorization(I);
    J:=Ideal(E1,ZBasis(I));
    assert Seqset(PrimesAbove(J)) eq { p[1] : p in fac };
    printf ".";

    // testing PlacesAboveRationalPrime: the primes above different rational primes are stored
    // separately
    _<x>:=PolynomialRing(Integers());
    A:=EtaleAlgebra((x^2+5)*(x^2+7));
    O:=MaximalOrder(A);
    P2:=PlacesAboveRationalPrime(A,2);
    P3:=PlacesAboveRationalPrime(A,3);
    assert Seqset(P2) eq Seqset(PrimesAbove(2*O));
    assert Seqset(P3) eq Seqset(PrimesAbove(3*O));
    assert Seqset(P2) ne Seqset(P3);
    assert PlacesAboveRationalPrime(A,2) eq P2; // the stored value
    assert Keys(A`PlacesAboveRationalPrime) eq {2,3};
    printf ".";

    // PrimeBelow: the prime S meet P of S below a prime P of an order T containing S. It is stored:
    // the same ideal object is returned again, also for an equal prime which is a different object.
    // Primes of the maximal order and of an intermediate order T, in an algebra with two
    // components.
    A:=EtaleAlgebra((x^2+5)*(x^3-2000));
    E:=EquationOrder(A);
    O:=MaximalOrder(A);
    T:=[ S : S in OverOrders(E) | S ne E and S ne O ][1];
    for P in PrimesAbove(Conductor(E)) do
        for U in [ O, T ] do
            for Q in PrimesAbove(U!!P) do
                Q2:=PrimeBelow(E,Q);
                assert Order(Q2) eq E and IsPrime(Q2) and Q2 eq OneIdeal(E) meet E!!Q;
                assert Q2 eq P; // E has a unique prime below Q, and P is a prime of E below Q
                assert IsIdentical(PrimeBelow(E,Q),Q2);
                Q3:=Ideal(U,ZBasis(Q)); // equal to Q, a different object
                assert not IsIdentical(Q,Q3) and IsIdentical(PrimeBelow(E,Q3),Q2);
            end for;
        end for;
    end for;
    // one ideal for each prime of E: the primes of E below the different primes of O and of T above
    // the same prime P of E are all the same object, which is P itself (PrimesAbove uses
    // PrimeBelow)
    for P in PrimesAbove(Conductor(E)) do
        below:=[ PrimeBelow(E,Q) : Q in PrimesAbove(O!!P) cat PrimesAbove(T!!P) ];
        assert forall{ Q2 : Q2 in below | IsIdentical(Q2,P) };
    end for;
    // P of T itself
    assert IsIdentical(PrimeBelow(T,PrimesAbove(Conductor(T))[1]),PrimesAbove(Conductor(T))[1]);
    try
        _:=PrimeBelow(O,PrimesAbove(Conductor(E))[1]); // O is not contained in E
        ok:=false;
    catch err
        ok:=true;
    end try;
    assert ok;
    printf ".";
    // at assertion level 2, the check of a factorization formed the prime powers with '^', which for
    // this ideal led to an infinite recursion (TwoGeneratingSet, CoprimeRepresentative, PrimesAbove,
    // Factorization) and a crash of Magma; found by the timing test PicardGroup_1.8 at level 2
    SetAssertions(2);
    _<x>:=PolynomialRing(Integers());
    A:=EtaleAlgebra(x^6-3*x^5-3*x^4+65*x^3-48*x^2-768*x+4096);
    assert [ DefiningPolynomial(K) : K in Components(A) ] eq [ x^2-7*x+16, x^4+4*x^3+9*x^2+64*x+256 ];
    zbS:=[ <[ 0, 0 ], [ -1/9, 1/48, 1/12, 11/144 ]>, <[ 0, 0 ], [ 0, -1/24, 1/6, -1/24 ]>,
           <[ 0, 0 ], [ 1/9, 3/16, 1/12, 19/144 ]>, <[ 0, 0 ], [ 1/9, -1/6, 0, 5/18 ]>,
           <[ 1, 0 ], [ 0, 0, 0, 0 ]>, <[ 0, 1 ], [ 0, 0, 0, 0 ]> ];
    gI:=[ <[ 2, 0 ], [ 0, 0, 0, 0 ]>, <[ 3, 1 ], [ 0, 0, 0, 0 ]>, <[ 0, 0 ], [ 1, 0, 0, 0 ]> ];
    S:=Order([ A!t : t in zbS ]);
    I:=Ideal(S,[ A!t : t in gI ]);
    TwoGeneratingSet(I);
    assert #Generators(I) le 2 and I eq Ideal(S,Generators(I));
    printf ".";
    SetAssertions(1);
    printf " all good!";

*/
