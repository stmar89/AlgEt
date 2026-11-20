/*
    TotRealPos_9.6: UnitGroup, TotallyRealUnitGroup and TotallyRealPositiveUnitGroup of all the overorders of nine
    orders in CM algebras of dimension 4 to 10: 1827 orders (OverOrders not timed).
    Why: the totally real subalgebra has to be computed, its unit group intersected with the unit group of each
    overorder through the unit group of the maximal order, and the totally positive units are found by solving a
    system over F_2 given by the signs at the real places. The orders: Z[F,q/F] for the Weil polynomials
    x^8-2x^7+7x^6-14x^5+40x^4-56x^3+112x^2-128x+256 (q = 2, see WKICM_4.7), x^8-x^7-3x^5-2x^4-15x^3-125x+625 (q = 5),
    x^10-x^8+2x^6+4x^4-8x^2+32 and x^10-x^9+4x^8-6x^7+8x^6-16x^5+16x^4-24x^3+32x^2-16x+32 (q = 2, products of CM
    fields) and x^6-3x^5+30x^4-175x^3+750x^2-1875x+15625 (q = 25); the equation orders Z[8 zeta_8] (x^4+4096),
    Z[2 zeta_7], Z[2 zeta_9] and the equation order of x^8+81 (q = 3). Replaces TotRealPos_9.3 (the 19 overorders of
    Z[F,2/F] for x^8+16, about 1 s with 9362a5e).
*/
AttachSpec(spec);
SetSeed(1);
SetClassGroupBounds("GRH");
_<x>:=PolynomialRing(Integers());
orders:=[];
for f in [ x^8-2*x^7+7*x^6-14*x^5+40*x^4-56*x^3+112*x^2-128*x+256, x^8-x^7-3*x^5-2*x^4-15*x^3-125*x+625,
           x^10-x^8+2*x^6+4*x^4-8*x^2+32, x^10-x^9+4*x^8-6*x^7+8*x^6-16*x^5+16*x^4-24*x^3+32*x^2-16*x+32,
           x^6-3*x^5+30*x^4-175*x^3+750*x^2-1875*x+15625 ] do
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    q:=Integers()!Round(ConstantCoefficient(f)^(2/Degree(f)));
    Append(~orders,Order([F,q/F]));
end for;
for f in [ x^4+4096, x^6+2*x^5+4*x^4+8*x^3+16*x^2+32*x+64, x^6+8*x^3+64, x^8+81 ] do
    Append(~orders,EquationOrder(EtaleAlgebra(f)));
end for;
oos:=[ OverOrders(R) : R in orders ];
assert [ #oo : oo in oos ] eq [ 420, 424, 176, 114, 126, 249, 119, 119, 80 ];
ranks:=[];
t0:=Cputime();
    for oo in oos do
        for S in oo do
            U,u:=UnitGroup(S);
            T:=TotallyRealUnitGroup(S);
            Tp:=TotallyRealPositiveUnitGroup(S);
            assert forall{ i : i in [1..Ngens(Tp)] | IsTotallyRealPositive(u(U!Tp.i)) };
            Append(~ranks,<AbelianInvariants(T),AbelianInvariants(Tp)>);
        end for;
    end for;
t:=Cputime(t0);
// the abelian invariants of the totally real units and of the totally positive ones
assert {* r : r in ranks *} eq {* <[ 2, 0 ], [ 0 ]>^^885, <[ 2, 0, 0 ], [ 0, 0 ]>^^238,
                                <[ 2, 0, 0, 0 ], [ 0, 0, 0 ]>^^80, <[ 2, 2, 0 ], [ 0 ]>^^436,
                                <[ 2, 2, 2, 0 ], [ 0 ]>^^161, <[ 2, 2, 2, 2, 0 ], [ 0 ]>^^27 *};
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
