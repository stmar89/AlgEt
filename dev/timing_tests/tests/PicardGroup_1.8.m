/*
    PicardGroup_1.8: PicardGroup (GRH) of all the overorders of six orders that are not equation orders, or whose
    overorders are not: 1482 orders in total (OverOrders not timed).
    Why: the whole PicardGroup machinery (residue rings, CRT, unit groups, gluing with Cl(O)) on many orders given by
    Z-bases, with different conductors and Picard groups, in six algebras: the order T given by a Z-basis in the
    algebra of the Weil polynomial x^6-3x^5-3x^4+65x^3-48x^2-768x+4096 (q = 16; an imaginary quadratic times a
    quartic CM field, Cl(O) = Z/2 x Z/4) and the order Z[F,16/F] in the same algebra; the equation orders of
    x^4+3^5 5^3 7^2 (CM quartic, Cl(O) = Z/2 x Z/4, 624 overorders), (x^2+5)(x^2+7)(x^2+11) and
    (x^2+1)(x^2+2)(x^2+3)(x^2+5) (products of imaginary quadratic fields); Z[F,2/F] for the Weil polynomial
    x^10-x^9+4x^8-6x^7+8x^6-16x^5+16x^4-24x^3+32x^2-16x+32 (q = 2, four components). Replaces PicardGroup_1.2
    (PicardGroup of T alone, 12 times, about 4 s with 9362a5e).
*/
AttachSpec(spec);
SetSeed(1);
SetClassGroupBounds("GRH");
_<x>:=PolynomialRing(Integers());
f16:=x^6-3*x^5-3*x^4+65*x^3-48*x^2-768*x+4096;
dataT:=[
    <[ 1, 0 ], [ 1/9, 5/6, 1, 41/18 ]>,
    <[ 0, 1 ], [ 0, 1, 0, 0 ]>,
    <[ 0, 0 ], [ 8/9, 11/6, 4/3, 133/18 ]>,
    <[ 0, 0 ], [ 0, 8/3, 7/3, 29/3 ]>,
    <[ 0, 0 ], [ 0, 0, 3, 3 ]>,
    <[ 0, 0 ], [ 0, 0, 0, 18 ]>
];
// the orders: Z[F,q/F] for f a Weil polynomial (q = constant term^(2/deg f)), the equation order, or T
orders:=[];
A:=EtaleAlgebra(f16);
assert [ DefiningPolynomial(K) : K in Components(A) ] eq [ x^2-7*x+16, x^4+4*x^3+9*x^2+64*x+256 ];
Append(~orders,Order([ A ! g : g in dataT ]));
A:=EtaleAlgebra(f16); F:=PrimitiveElement(A);
Append(~orders,Order([F,16/F]));
for f in [ x^4+3^5*5^3*7^2, (x^2+5)*(x^2+7)*(x^2+11), (x^2+1)*(x^2+2)*(x^2+3)*(x^2+5) ] do
    Append(~orders,EquationOrder(EtaleAlgebra(f)));
end for;
A:=EtaleAlgebra(x^10-x^9+4*x^8-6*x^7+8*x^6-16*x^5+16*x^4-24*x^3+32*x^2-16*x+32); F:=PrimitiveElement(A);
Append(~orders,Order([F,2/F]));
oos:=[ OverOrders(R) : R in orders ];
assert [ #oo : oo in oos ] eq [ 208, 144, 624, 276, 116, 114 ];
nums:=[];
t0:=Cputime();
    for oo in oos do
        Append(~nums,[ #PicardGroup(S : GRH:=true) : S in oo ]);
    end for;
t:=Cputime(t0);
// the sums of the orders of the Picard groups of the overorders
assert [ &+n : n in nums ] eq [ 8256, 38176, 984096, 5236, 520, 437 ];
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
