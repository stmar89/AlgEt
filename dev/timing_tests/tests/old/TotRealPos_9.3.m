/*
    TotRealPos_9.3: TotallyRealUnitGroup and TotallyRealPositiveUnitGroup of the 19 overorders of R = Z[F,2/F] in the
    CM octic field of f = x^8+16 (OverOrders not timed).
    Why: the CM field Q(zeta_16)'s cousin of x^8+16 has a totally real subfield of degree 4 whose unit group (rank 3)
    has to be computed and intersected with the unit group of each overorder, through the unit group of the maximal
    order (rank 3, with roots of unity of order 16); the totally positive units are found by solving a system over
    F_2 given by the signs at the 4 real places. The same overorders as in WKICM_4.1.
*/
AttachSpec(spec);
SetSeed(1);
SetClassGroupBounds("GRH");
_<x>:=PolynomialRing(Integers());
f:=x^8+16;
A:=EtaleAlgebra(f);
F:=PrimitiveElement(A);
R:=Order([F,2/F]);
oo:=OverOrders(R);
assert #oo eq 19;
t0:=Cputime();
    for S in oo do
        U,u:=UnitGroup(S);
        T:=TotallyRealUnitGroup(S);
        Tp:=TotallyRealPositiveUnitGroup(S);
        assert forall{ i : i in [1..Ngens(Tp)] | IsTotallyRealPositive(u(U!Tp.i)) };
    end for;
t:=Cputime(t0);
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
