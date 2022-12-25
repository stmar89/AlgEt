/* vim: set syntax=magma :*/
    "-------------------------------------------------------------";
    "-------------------------------------------------------------";
    "Slow tests, maily to compare timings (on diophantus)";
    "Code with profiler is commented out";
    "If the timing is off by more than 10%% with respect to the previous best (achieved on diophantus) then a warning is printed.";
    time_start_slow_tests:=Cputime();

    "-------------------------------------------------------------";
    "-------------------------------------------------------------";
    "------------------------Picard Groups------------------------";
    "-------------------------------------------------------------";
    "-------------------------------------------------------------";
    "Picard group for AlgEtQOrd used to have very incosistent timings because of how certain randomized choices in CoprimeRepresentatives. Now the timings are consistent and better than RngOrd.";
    //picard groups
	AttachSpec("~/packages_github/AlgEt/spec");
    SetClassGroupBounds("GRH");
	_<x>:=PolynomialRing(Integers());
    f:=x^4-1000*x^3-1000*x^2-1000*x-1000;
    t_tot:=0;
    "RngOrd";
    K2:=NumberField(f);
    E2:=EquationOrder(K2);
    time P2,p2:=PicardGroup(E2);
    for i in [1..5] do
        printf "%o-th round.\n",i;
        "AlgEt";
        K1:=EtaleAlgebra(f);
        E1:=EquationOrder(K1);
        // no profiler
        t0:=Cputime();
            time P1,p1:=PicardGroup(E1 : GRH:=true );
        t1:=Cputime(t0);
        t_tot+:=t1;
        assert #P1 eq #P2;
    end for;
    printf "\n";
    printf "Current running time = %o \n",t_tot;
    tprevbest:=109.890; //on diophantus
    if Abs(t_tot - tprevbest) gt 0.1*tprevbest then
        if t_tot lt tprevbest then
            printf "PicardGroup for %o got faster. Update the previous best known time\n",f;
        elif t_tot gt tprevbest then
            printf "PicardGroup for %o got slower. The previous code was better.\n",f;
        end if;
    end if;


    // with profiler
    // K1:=EtaleAlgebra(f);
    // E1:=EquationOrder(K1);
    // SetProfile(true);
    //     P1,p1:=PicardGroup(E1);
    // SetProfile(false);
    // ProfilePrintByTotalTime(ProfileGraph());

    "-------------------------------------------------------------";
    "The following order used to trigger a bug in CRT. now fixed";
    clear;
    AttachSpec("~/packages_github/AlgEt/spec");
    _<x>:=PolynomialRing(Integers());
    f:=x^6 - 3*x^5 - 3*x^4 + 65*x^3 - 48*x^2 - 768*x + 4096;
    A:=EtaleAlgebra(f);
    gensT:=[
        <[ 1, 0 ], [ 1/9, 5/6, 1, 41/18 ]>,
        <[ 0, 1 ], [ 0, 1, 0, 0 ]>,
        <[ 0, 0 ], [ 8/9, 11/6, 4/3, 133/18 ]>,
        <[ 0, 0 ], [ 0, 8/3, 7/3, 29/3 ]>,
        <[ 0, 0 ], [ 0, 0, 3, 3 ]>,
        <[ 0, 0 ], [ 0, 0, 0, 18 ]>
    ];
    gensT:=[ A ! g : g in gensT ];
    T:=Order(gensT);
    t0:=Cputime();
        _:=#PicardGroup(T);
    t1:=Cputime(t0);
    printf "Current running time = %o \n",t1;
    tprevbest:=13.7; //on diophantus
    if Abs(t1 - tprevbest) gt 0.1*tprevbest then
        if t1 lt tprevbest then
            printf "PicardGroup for %o got faster. Update the previous best known time\n",f;
        elif t1 gt tprevbest then
            printf "PicardGroup for %o got slower. The previous code was better.\n",f;
        end if;
    end if;

	
    "-------------------------------------------------------------";
    "-------------------------------------------------------------";
    "--------------------------OverOrders-------------------------";
    "-------------------------------------------------------------";
    "-------------------------------------------------------------";
    "A quick one.";
    clear;
	_<x>:=PolynomialRing(Integers());
    f:=x^4-10000*x^3-10000*x^2-10000*x-10000; 
    AttachSpec("~/packages_github/AlgEt/spec");
    // no profiler
    Aet:=EtaleAlgebra(f);
    Eet:=EquationOrder(Aet);
    t0:=Cputime();
        ooet:=FindOverOrders(Eet);
    t1:=Cputime(t0); t1;
    tprevbest:=13.7; //on diophantus
    if Abs(t1 - tprevbest) gt 0.1*tprevbest then
        if t1 lt tprevbest then
            printf "FindOverOrders for %o got faster. Update the previous best known time\n",f;
        elif t1 gt tprevbest then
            printf "FindOverOrders for %o got slower. The previous code was better.\n",f;
        end if;
    end if;
    assert #ooet eq 297;

    // Aet:=EtaleAlgebra(f);
    // Eet:=EquationOrder(Aet);
    // SetProfile(true);
    // ooet:=FindOverOrders(Eet);
    // SetProfile(false);
    // G2:=ProfileGraph();
    // ProfilePrintByTotalTime(ProfilePruneGraphByTime(G2,30));

    "-------------------------------------------------------------";
    "A bigger example.";
    clear;
    AttachSpec("~/packages_github/AlgEt/spec");
    // no profiler
	_<x>:=PolynomialRing(Integers());
    f:=x^6 - 3*x^5 - 3*x^4 + 65*x^3 - 48*x^2 - 768*x + 4096;
    Aet:=EtaleAlgebra(f);
    Eet:=EquationOrder(Aet);
    t0:=Cputime();
        ooet:=FindOverOrders(Eet);
    t1:=Cputime(t0);
    printf "Current running time = %o \n",t1;
    tprevbest:=693; // on diophantus
    if Abs(t1 - tprevbest) gt 0.1*tprevbest then
        if t1 lt tprevbest then
            printf "FindOverOrders for %o got faster. Update the previous best known time\n",f;
        elif t1 gt tprevbest then
            printf "FindOverOrders for %o got slower. The previous code was better.\n",f;
        end if;
    end if;
    assert #ooet eq 3312; 

    /*
    // with profiler
    Aet:=EtaleAlgebra(f);
    Eet:=EquationOrder(Aet);
    SetProfile(true);
    ooet:=FindOverOrders(Eet);
    SetProfile(false);
    G2:=ProfileGraph();
    ProfilePrintByTotalTime(ProfilePruneGraphByTime(G2,30));
    */

    "-------------------------------------------------------------";
    "-------------------------------------------------------------";
    "------------------------WkEquivClasses-----------------------";
    "-------------------------------------------------------------";
    "-------------------------------------------------------------";

    "A very quick example.";
    clear;
	_<x>:=PolynomialRing(Integers());
    f:=x^8+16; 
    AttachSpec("~/packages_github/AlgEt/spec");

    // no profiler
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    R:=Order([F,2/F]); 
    oo:=FindOverOrders(R : populateoo_in_oo:=true); //to avoid any computation related to the OverOrders
    t0:=Cputime();
        #WKICM(R);
    t1:=Cputime(t0);
    printf "Current running time = %o \n",t1;
    tprevbest:=4.5; // diophantus
    if Abs(t1 - tprevbest) gt 0.1*tprevbest then
        if t1 lt tprevbest then
            printf "WKICM for %o got faster. Update the previous best known time\n",f;
        elif t1 gt tprevbest then
            printf "WKICM for %o got slower. The previous code was better.\n",f;
        end if;
    end if;
    /* 
    // with profiler
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    R:=Order([F,2/F]); 
    oo:=FindOverOrders(R);
    SetProfile(true);
        #WKICM(R);
    SetProfile(false);
    G2:=ProfileGraph();
    ProfilePrintByTotalTime(ProfilePruneGraphByTime(G2,30));
    */
    "-------------------------------------------------------------";
    "A much bigger example: the interesting overorders, with CohenMacaulyType > 2, have been already computed and saved in a special file. We compute the WkICM_bar for this order, and compare timings. The orders are sorted from the fastest to the slowest.";
    clear;
	_<x>:=PolynomialRing(Integers());
    f:=x^8+16; 
    AttachSpec("~/packages_github/AlgEt/spec");
    A:=EtaleAlgebra(f);
    R:=EquationOrder(A);
    ooR:=FindOverOrders(R);
    data:=eval(Read("~/packages_github/AlgEt/dev/input_big_test_WKICM.txt"));
    // data contains a Seq with entries <S,T,ff,out> where
    // S is the ZBasis of the order S
    // T (resp. ff) is the ZBasis of the overorder T of S (resp ff=(S:T))
    // out is the number of intermedite S-modules between T and ff.
    // data is already sorted wrt to the size of out.
    oo:={@ @};
    for x in data do
        S:=Order([ A!z:z in x[1]]);
        S`OverOrders:={@ T :T in ooR | S subset T @};
        Include(~oo,S);
    end for;
    // "IntermediateIdeals";
    // SetProfile(true);
    out_new:=[]; //oo[15] takes ~100 seconds
    for i->S in oo do
        t0:=Cputime();
        delete S`WKICM_bar;
        N:=#WKICM_bar(S : Method:="IntermediateIdeals");
        t1:=Cputime(t0);
        out:=<i,N,t1>; out;
        Append(~out_new,out);
    end for;
    // out on diophantus
    out_prev:=[
    <1, 4, 3.130>,
    <2, 4, 2.300>,
    <3, 2, 2.180>,
    <4, 4, 2.810>,
    <5, 2, 1.840>,
    <6, 4, 2.890>,
    <7, 4, 2.930>,
    <8, 4, 2.540>,
    <9, 2, 1.770>,
    <10, 4, 2.520>,
    <11, 4, 2.140>,
    <12, 4, 2.450>,
    <13, 4, 2.190>,
    <14, 6, 24.570>,
    <15, 6, 92.750>,
    <16, 6, 96.850>,
    <17, 6, 95.610>,
    <18, 6, 216.540>,
    <19, 6, 772.750> 
    ];
    // SetProfile(false);
    // ProfilePrintByTotalTime(ProfileGraph() : Max:=30);
    assert forall{ i : i in [1..#out_prev] | out_prev[i,2] eq out_new[i,2] };
    tprevbest:=&+[ o[3] : o in out_prev ];
    t1:=&+[ o[3] : o in out_new ];
    if Abs(t1 - tprevbest) gt 0.1*tprevbest then
        if t1 lt tprevbest then
            printf "WKICM for %o got faster. Update the previous best known time\n",f;
        elif t1 gt tprevbest then
            printf "WKICM for %o got slower. The previous code was better.\n",f;
        end if;
    end if;
    // "LowIndexProcess";
    // outputLowIndexProcess:=[];
    // for i->S in oo[1..14] do //the last five are very slow.
    //     t0:=Cputime();
    //     N:=#WKICM_bar(S : Method:="LowIndexProcess");
    //     t1:=Cputime(t0);
    //     out:=<i,N,t1>; out;
    //     Append(~outputLowIndexProcess,out);
    //     assert N eq outputIntermediateIdeals[i][2];
    // end for;
    "-------------------------------------------------------------";
    printf "The whole set of slow tests required %o\n",Cputime(time_start_slow_tests);
