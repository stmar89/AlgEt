/*
    IsPrincipal_6.2: IsPrincipal (GRH) on two principal ideals of the equation order E of f = x^8 + 3^8 7, each with
    the random seeds 1, 2 and 3 on a fresh copy (PicardGroup not timed).
    Why: the order of IsPrincipal_6.1. The generator of a principal ideal I is obtained from a generator of IO and an
    arbitrary coset representative in the unit group of O, with exponents that can be in the hundreds of thousands,
    and the running time depends on the random coprime representative chosen inside IsPrincipal (hence the three
    seeds): from a second to more than 40 minutes per ideal with a2f3804. The benchmark for plan item B1 (cache the
    unit data and reduce the exponents modulo the units of E), which should make this test take seconds; until then
    it times out (40 minutes) with both v1.0.2 and a2f3804.
    Inputs (fixed, see README): two ideals p(g)*p(-g), for random classes g, whose generator the short-element
    shortcut does not find, computed with a2f3804 and Magma 2.30-1 on 2026-11-10 and stored as Z-bases.
*/
AttachSpec(spec);
SetSeed(1);
SetClassGroupBounds("GRH");
_<x>:=PolynomialRing(Integers());
// Fixed inputs (see README): a checksum of the coordinates of a sequence of elements, printed and asserted below.
checksum:=function(elts)
    s:=0;
    for z in elts do
        for c in Components(z) do
            for q in Eltseq(c) do
                s:=(1000003*s+7*Numerator(q)+Denominator(q)) mod (2^61-1);
            end for;
        end for;
    end for;
    return s;
end function;
// fixed inputs (one Z-basis per ideal), printed with PrintSeqAlgEtQElt
data_principal:=[
    [<[194157,568548,104703,118299,-617171,-323040,-158365,-305695]>,<[57591,-100818,502461,58662,121677,663114,
    451751,321890]>,<[-122472,-469152,670743,129111,-317238,-448010,410547,-143961]>,<[-364257,938772,304476,194568,
    182005,-635854,534752,615926]>,<[-486000,-68760,-161868,-534903,-842582,738420,-571202,139115]>,<[-930447,
    -529155,199839,557604,-217255,610523,10003,-252316]>,<[-124416,339831,-748479,-163671,205612,-281025,1041299,
    -44163]>,<[-613818,880821,428580,-251169,579834,8411,264588,-811343]>],
    [<[-291616,-527476,-294127,-601997,-289030,141776,-465859,221541]>,<[-461586,231932,107251,-400,-206424,-324548,
    992981,52954]>,<[-363686,518895,-391236,-655899,-407947,76938,177447,-282450]>,<[-123184,-510562,333319,211392,
    -637354,345790,487653,32618]>,<[727984,-117056,28979,181208,-510348,-157094,467663,738780]>,<[-28740,-215835,
    -343052,-366689,373226,1065409,213302,-292021]>,<[175063,564070,874949,-135820,70289,445738,89587,487912]>,
    <[-506493,132049,68095,650421,-93683,689447,-596407,322747]>]
];
f:=x^8+3^8*7;
A:=EtaleAlgebra(f);
E:=EquationOrder(A);
copies:=[ [ Ideal(E,[ A!t : t in zb ]) : seed in [1..3] ] : zb in data_principal ];
chk:=checksum([ A!t : t in &cat data_principal ]);
printf "### INPUT %o\n",chk;
assert chk eq 47279510205790644;
P,p:=PicardGroup(E : GRH:=true);
assert #P eq 1259712;
t0:=Cputime();
    for Is in copies do
        for seed->I in Is do
            SetSeed(seed);
            test,gen:=IsPrincipal(I : GRH:=true);
            assert test and gen*E eq I;
        end for;
    end for;
t:=Cputime(t0);
printf "### TIMING %o %o\n",name,t;
print "### DONE";
quit;
