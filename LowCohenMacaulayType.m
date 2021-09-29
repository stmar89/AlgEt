/* vim: set syntax=magma :*/

//freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////


/*TODO:

*/

declare attributes AlgEtOrd:NonGorensteinPrimes;

//------------
// NonGorensteinPrimes
//------------

intrinsic NonGorensteinPrimes(S::AlgEtOrd)->SeqEnum,SeqEnum
{ Given an order S it returns two sequences the first containing the primes at which S is locally not Gorenstein and the second containing the CohenMacaulay types of S at this primes, that is, the dimension of S^t/PS^t over S/P, where S^t is the TraceDualIdeal of S}
    if not assigned S`NonGorensteinPrimes then
        if IsGorenstein(S) then
            S`NonGorensteinPrimes:=<[],[]>
        else
            St:=TraceDualIdeal(S);
            pp0:=PrimesAbove(Conductor(S));
            pp:=[];
            dPs:=[];
            for P in pp0 do
                k:=Integers() ! Index(S,P);
                v:=Integers() ! Index(St,P*St);
                dP:=Ilog(k,v);
                if dP gt 1 then // don't want P such that S_P is Gorenstein
                    Append(~pp,P);
                    Append(~dPs,dP);
                end if;
            end for;
            assert Min(dPs) ne 1;
            S`NonGorensteinPrimes:=< pp,dPs >;
        end if;
    end if;
    return Explode(S`NonGorensteinPrimes);
end intrinsic;



/////////////////////////////////////////////////////
// functions for orders that locally have CM-type \leq 2
/////////////////////////////////////////////////////
// In a so-far-unpublished note I have proved that 
// if S^t/P*S^t has dimension 2 over S/P then
// locally at P the only fractional ideals with 
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
        pows_hat:=[ &*[ pows[j] : j in [1..#pp] | j ne i ] : i in [1..#pp] ];
        Is:=[ seq : seq in CartesianProduct([[OneIdeal(S),I] : i in [1..#pp]]) ];
        out:=[];
        for I in Is do
            L:=&+[ (I[i]+pows[i]) * pows_hat[i] : i in [1..#pows_hat]];
            Append(~out,L);
        end for;
        return out;
     end if;
end function;



/* TEST

*/
