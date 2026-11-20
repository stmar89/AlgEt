/*
    ICM_5.8: ICM (GRH) of S = Z + 2O in 14 products of three imaginary quadratic fields Q(sqrt(-d1)) x Q(sqrt(-d2)) x
    Q(sqrt(-d3)): 46460 ideal classes in total.
    Why: 2 is ramified, split or inert in each factor depending on d mod 8, so O/2O has dimension 6 with different
    structures and the unique singular prime P = 2O of S has S/P = F_2; S has Cohen-Macaulay type 5 and from 31 to
    104 overorders, WKICM(S) has from 126 to 823 classes (the general branch), and Cl(O) goes from trivial to
    Z/6 (for d = 2, 5, 23). The test runs the whole chain WKICM, PicardGroup of the
    overorders, ICM on 14 different orders. Replaces ICM_5.3 (only d = 5, 7, 11, 3712 classes, about 12 s with
    9362a5e).
*/
AttachSpec(spec);
SetSeed(1);
SetClassGroupBounds("GRH");
_<x>:=PolynomialRing(Integers());
// < [d1,d2,d3], number of ideal classes of Z + 2O >
data:=[ <[5,7,11],3712>, <[1,2,3],408>, <[1,2,5],2104>, <[3,5,7],1280>, <[2,7,11],1856>, <[1,3,7],357>,
        <[5,11,13],8416>, <[2,3,5],1456>, <[7,11,15],3292>, <[1,5,13],4208>, <[3,11,19],918>, <[2,5,23],10032>,
        <[5,7,15],5952>, <[1,7,23],2469> ];
assert &+[ d[2] : d in data ] eq 46460;
t:=0;
for d in data do
    ds,n:=Explode(d);
    A:=EtaleAlgebra(&*[ x^2+e : e in ds ]);
    O:=MaximalOrder(A);
    S:=Order([One(A)] cat [2*z : z in ZBasis(O)]);
    t0:=Cputime();
        icm:=ICM(S : GRH:=true);
    t+:=Cputime(t0);
    assert #icm eq n;
end for;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
