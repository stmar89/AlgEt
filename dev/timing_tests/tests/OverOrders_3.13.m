/*
    OverOrders_3.13: OverOrders of S = Z + pO for 7 algebras where O/pO is far from a field: 989 overorders in total.
    Why: S has a unique singular prime P = pO, with S/P = F_p and (P:P)/P = O/pO of dimension n-1, so S has
    Cohen-Macaulay type n-1, the maximum possible, and its local lattice of overorders is deep. The algebras:
    Q(2^(1/9)) (2 totally ramified, e = 9; the slowest), the product of Q(i), Q(sqrt(-2)), Q(sqrt(-5)) and
    Q(sqrt(-7)) (2 ramified in three components, split in the last), the totally real fields Q(zeta_32)^+ (2 totally
    ramified) and Q(zeta_17)^+ (2 splits into two primes of residue degree 4), Q(3^(1/6)) and Q(zeta_9) (3 totally
    ramified) and Q(5^(1/5)) (5 totally ramified). The type-2 enumeration of subspaces of T/P in
    MinimalOverOrdersAtPrime is the cost. Replaces OverOrders_3.8 (Z + 2O in Q(zeta_16), about 5 s with 9362a5e).
*/
AttachSpec(spec);
SetSeed(1);
_<x>:=PolynomialRing(Integers());
data:=[ <x^9-2,2,193>, <(x^2+1)*(x^2+2)*(x^2+5)*(x^2+7),2,516>, <x^8-8*x^6+20*x^4-16*x^2+2,2,110>,
        <x^8+x^7-7*x^6-6*x^5+15*x^4+10*x^3-10*x^2-4*x+1,2,16>, <x^6+3,3,42>, <x^6+x^3+1,3,42>, <x^5-5,5,15> ];
t:=0;
for d in data do
    f,p,n:=Explode(d);
    A:=EtaleAlgebra(f);
    O:=MaximalOrder(A);
    S:=Order([One(A)] cat [p*z : z in ZBasis(O)]);
    t0:=Cputime();
        oo:=OverOrders(S);
    t+:=Cputime(t0);
    assert #oo eq n;
end for;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
