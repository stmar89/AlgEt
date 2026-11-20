/*
    PicardGroup_1.7: PicardGroup (GRH) of the equation orders E of f = x^4 - 1000k(x^3+x^2+x+1) for k = 1..8, 10, 12,
    each on 3 fresh copies of the algebra, with the seeds 1, 2 and 3 (the result is cached in the order).
    Why: quartic fields of signature (2,1) with big class groups (Z/3548 for k = 1, Z/4 x Z/14796 for k = 3, Z/9764
    for k = 10); [O:E] is divisible by 2^3 5^3 and E has singular primes above 2, 5 and the primes dividing k, so
    the discrete logarithms in Pic(O) and the gluing with (O/f)^* are not trivial, and #Pic(E) goes from 3548000 to
    3801088000. The running time depends a lot on the random choices in CoprimeRepresentative (from 12 s to 49 s for
    k = 3, depending on the seed), hence the three seeds. The values of #Pic(E) agree with Magma's PicardGroup of the
    RngOrd E (checked in November 2026, not in the test since it takes up to 2 minutes per order). Replaces
    PicardGroup_1.1 (k = 1 on 5 fresh copies, about 2 s with 9362a5e).
*/
AttachSpec(spec);
SetSeed(1);
SetClassGroupBounds("GRH");
_<x>:=PolynomialRing(Integers());
// < k, #Pic(E) >
data:=[ <1,3548000>, <2,10240000>, <3,59184000>, <4,119808000>, <5,191724000>, <6,329088000>, <7,484096000>,
        <8,545792000>, <10,1249792000>, <12,3801088000> ];
t:=0;
for d in data do
    k,n:=Explode(d);
    f:=x^4-1000*k*(x^3+x^2+x+1);
    for seed in [1,2,3] do
        SetSeed(seed);
        A:=EtaleAlgebra(f);
        E:=EquationOrder(A);
        t0:=Cputime();
            P,p:=PicardGroup(E : GRH:=true);
        t+:=Cputime(t0);
        assert #P eq n;
    end for;
end for;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
