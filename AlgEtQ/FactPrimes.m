/* vim: set syntax=magma :*/

freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose FactPrimes, 3;

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
{Check if the order is Bass at the prime ideal P.}
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
{Check if the order is Bass.}
// we compute the maximal order O and check if O/PO is at most 2-dimensional over S/P for every singular prime P
// This coincides with the usual definition since O has the maximal number of generators as a fractional S ideal and S is Bass iff every ideal can be generated by at most 2-elements.
// see Hofman,Sircana, On the computations of over-orders, Definition 5.23
    if not assigned S`IsBass then
        if IsMaximal(S) then 
            S`IsBass:=true;
        else
            O:=MaximalOrder(Algebra(S));
            ff:=Conductor(S);
            sing:=PrimesAbove(ff);
            S`IsBass:=forall{P:P in sing | IsBassAtPrime(S,P)};
        end if;
    end if;
    return S`IsBass;
end intrinsic;

intrinsic IsGorensteinAtPrime(S::AlgEtQOrd,P::AlgEtQIdl) -> BoolElt
{Check if the order is Gorenstein at the prime ideal P.}
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


/* TESTS

    printf "### Testing Primes and Factorizaton:";
    Attach("~/packages_github/AlgEtQ/AlgEtQ.m");
    Attach("~/packages_github/AlgEtQ/Elt.m");
    Attach("~/packages_github/AlgEtQ/Ord.m");
    Attach("~/packages_github/AlgEtQ/TraceNorm.m");
    Attach("~/packages_github/AlgEtQ/Idl.m");
    Attach("~/packages_github/AlgEtQ/WkTesting.m");
    Attach("~/packages_github/AlgEtQ/FactPrimes.m");
    SetAssertions(2);

    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16);
    A:=EtaleAlgebra(f);
    E1:=EquationOrder(A);
    ff:=Conductor(E1);
    _:=PrimesAbove(Conductor(E1));
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
    facs:=[ Factorization(I) : I in ids | I ne OneIdeal(E1) ];
    printf ".";
    SetAssertions(1);
    printf " all good!\n"; 

*/
