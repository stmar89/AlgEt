/*
    PicardGroup_1.2: PicardGroup (GRH) of an order T given by a Z-basis (not monogenic) in the sextic algebra of
    f = x^6-3x^5-3x^4+65x^3-48x^2-768x+4096, a Weil polynomial for q = 16, computed on 12 fresh copies of the
    algebra (the result is cached in the order).
    Why: f = (x^2-7x+16)(x^4+4x^3+9x^2+64x+256): the algebra is the product of an imaginary quadratic field and a
    quartic CM field, with Cl(O) = Z/2 x Z/4. [O:T] = 55296 = 2^11 3^3; T has 3 singular primes (two above 2 with
    residue field F_2, one above 3), is Gorenstein and has 208 overorders; #Pic(T) = 192. T used to trigger a bug in
    ChineseRemainderTheorem, now fixed: a quick test of the whole PicardGroup machinery (residue rings, CRT, unit
    groups, gluing with Cl(O)) on an order that is not an equation order.
*/
AttachSpec(spec);
SetSeed(1);
SetClassGroupBounds("GRH");
_<x>:=PolynomialRing(Integers());
f:=x^6 - 3*x^5 - 3*x^4 + 65*x^3 - 48*x^2 - 768*x + 4096;
dataT:=[
    <[ 1, 0 ], [ 1/9, 5/6, 1, 41/18 ]>,
    <[ 0, 1 ], [ 0, 1, 0, 0 ]>,
    <[ 0, 0 ], [ 8/9, 11/6, 4/3, 133/18 ]>,
    <[ 0, 0 ], [ 0, 8/3, 7/3, 29/3 ]>,
    <[ 0, 0 ], [ 0, 0, 3, 3 ]>,
    <[ 0, 0 ], [ 0, 0, 0, 18 ]>
];
t:=0;
for i in [1..12] do
    A:=EtaleAlgebra(f);
    T:=Order([ A ! g : g in dataT ]);
    t0:=Cputime();
        P,p:=PicardGroup(T : GRH:=true);
    t+:=Cputime(t0);
    assert #P eq 192;
end for;
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
