/* vim: set syntax=magma :*/

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, stefano.marseglia89@gmail.com
// https://stmar89.github.io/index.html
// 
// Distributed under the terms of the GNU Lesser General Public License (L-GPL)
//      http://www.gnu.org/licenses/
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 3.0 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301  USA
// 
// Copyright 2023, S. Marseglia
/////////////////////////////////////////////////////

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
	AttachSpec("../spec");
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
    tprevbest:=47.880; //on diophantus
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
    
    AttachSpec("../spec");
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
    tprevbest:=0.5; //on diophantus
    if Abs(t1 - tprevbest) gt tprevbest then //so fast than the small random component in SmallRepresentative triggers the warning
        if t1 lt tprevbest then
            printf "PicardGroup for %o got faster. Update the previous best known time\n",f;
        elif t1 gt tprevbest then
            printf "PicardGroup for %o got slower. The previous code was better.\n",f;
        end if;
    end if;
	
    "-------------------------------------------------------------";
    "-------------------------------------------------------------";
    "-----------------OverOrders and GraphOverOrders--------------";
    "-------------------------------------------------------------";
    "-------------------------------------------------------------";
    "A quick one.";
    
    AttachSpec("../spec");
	_<x>:=PolynomialRing(Integers());
    f:=x^4-10000*x^3-10000*x^2-10000*x-10000; 
    // no profiler
    Aet:=EtaleAlgebra(f);
    Eet:=EquationOrder(Aet);
    t0:=Cputime();
        ooet:=FindOverOrders(Eet);
    t1:=Cputime(t0);
    printf "Current running time = %o \n",t1;
    tprevbest:=1.1; //on diophantus
    if Abs(t1 - tprevbest) gt 0.1*tprevbest then
        if t1 lt tprevbest then
            printf "FindOverOrders for %o got faster. Update the previous best known time\n",f;
        elif t1 gt tprevbest then
            printf "FindOverOrders for %o got slower. The previous code was better.\n",f;
        end if;
    end if;
    assert #ooet eq 297;

    t0:=Cputime();
        ooet:=GraphOverOrders(Eet);
    t1:=Cputime(t0);
    printf "Current running time = %o \n",t1;
    tprevbest:=0.3; //on diophantus
    if Abs(t1 - tprevbest) gt 1*tprevbest then
        if t1 lt tprevbest then
            printf "GraphOverOrders for %o got faster. Update the previous best known time\n",f;
        elif t1 gt tprevbest then
            printf "GraphOverOrders for %o got slower. The previous code was better.\n",f;
        end if;
    end if;

    // Aet:=EtaleAlgebra(f);
    // Eet:=EquationOrder(Aet);
    // SetProfile(true);
    // ooet:=FindOverOrders(Eet);
    // SetProfile(false);
    // G2:=ProfileGraph();
    // ProfilePrintByTotalTime(ProfilePruneGraphByTime(G2,30));

    "-------------------------------------------------------------";
    "A bigger example.";
    
    AttachSpec("../spec");
    // no profiler
	_<x>:=PolynomialRing(Integers());
    f:=x^6 - 3*x^5 - 3*x^4 + 65*x^3 - 48*x^2 - 768*x + 4096;
    Aet:=EtaleAlgebra(f);
    Eet:=EquationOrder(Aet);
    t0:=Cputime();
        ooet:=FindOverOrders(Eet);
    t1:=Cputime(t0);
    printf "Current running time = %o \n",t1;
    tprevbest:=15; // on diophantus
    if Abs(t1 - tprevbest) gt 0.1*tprevbest then
        if t1 lt tprevbest then
            printf "FindOverOrders for %o got faster. Update the previous best known time\n",f;
        elif t1 gt tprevbest then
            printf "FindOverOrders for %o got slower. The previous code was better.\n",f;
        end if;
    end if;
    assert #ooet eq 3312; 

    t0:=Cputime();
        ooet:=GraphOverOrders(Eet);
    t1:=Cputime(t0);
    printf "Current running time = %o \n",t1;
    tprevbest:=30.5; //on diophantus
    if Abs(t1 - tprevbest) gt 0.1*tprevbest then
        if t1 lt tprevbest then
            printf "GraphOverOrders for %o got faster. Update the previous best known time\n",f;
        elif t1 gt tprevbest then
            printf "GraphOverOrders for %o got slower. The previous code was better.\n",f;
        end if;
    end if;
    
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
    "A big example.";

    AttachSpec("../spec");
    P<x>:=PolynomialRing(Integers());
    f:=x^6 + 8*x^5 + 50*x^4 + 200*x^3 + 1250*x^2 + 5000*x + 15625;
    A:=EtaleAlgebra(f);
    Eet:=EquationOrder(A); 
    t0:=Cputime();
        ooet:=FindOverOrders(Eet);
    t1:=Cputime(t0);
    assert #ooet eq 3200;
    printf "Current running time = %o \n",t1;
    tprevbest:=21; // on diophantus
    if Abs(t1 - tprevbest) gt 0.1*tprevbest then
        if t1 lt tprevbest then
            printf "FindOverOrders for %o got faster. Update the previous best known time\n",f;
        elif t1 gt tprevbest then
            printf "FindOverOrders for %o got slower. The previous code was better.\n",f;
        end if;
    end if;
    t0:=Cputime();
        ooet:=GraphOverOrders(Eet);
    t1:=Cputime(t0);
    printf "Current running time = %o \n",t1;
    tprevbest:=37; //on diophantus
    if Abs(t1 - tprevbest) gt 0.1*tprevbest then
        if t1 lt tprevbest then
            printf "GraphOverOrders for %o got faster. Update the previous best known time\n",f;
        elif t1 gt tprevbest then
            printf "GraphOverOrders for %o got slower. The previous code was better.\n",f;
        end if;
    end if;
        
    "-------------------------------------------------------------";
    "-------------------------------------------------------------";
    "------------------------WkEquivClasses-----------------------";
    "-------------------------------------------------------------";
    "-------------------------------------------------------------";

    "A very quick example.";
    
	_<x>:=PolynomialRing(Integers());
    f:=x^8+16; 
    AttachSpec("../spec");

    // no profiler
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    R:=Order([F,2/F]); 
    oo:=FindOverOrders(R : populateoo_in_oo:=true); //to avoid any computation related to the OverOrders
    t0:=Cputime();
        _:=#WKICM(R);
    t1:=Cputime(t0);
    printf "Current running time = %o \n",t1;
    tprevbest:=2.1; // diophantus
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
    "An example in which several overorders have more than one singular prime";

    P<x>:=PolynomialRing(Integers());
    f:=x^6 + 8*x^5 + 50*x^4 + 200*x^3 + 1250*x^2 + 5000*x + 15625;
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    gensS:=[
    <[1,0],[1,0],[1,0]>,
    <[0,1],[0,1],[0,1]>,
    <[0,-2],[0,4],[0,-2]>,
    <[8,0],[-2,0],[-5,0]>,
    <[0,-6],[0,-2],[0,7]>,
    <[2,0],[-18,0],[15,0]>
    ];
    S:=Order([A!g: g in gensS]);
    t0:=Cputime();
        assert #WKICM_bar(S) eq 7;   
    t_curr:=Cputime(t0);
    t_prev_best:=11.9;
    "Current running time: ",t_curr;
    if Abs((t_curr - t_prev_best)/t_prev_best) gt 0.1 then 
        print "The current timing is different from the previous best known one. UPDATE!"; 
    end if;

    "-------------------------------------------------------------";
    "An example in which the order has a unique singular prime";

    // R has one singular prime
    P<x>:=PolynomialRing(Integers());
    f:=x^8+16;
    A:=EtaleAlgebra(f);
    R:=EquationOrder(A);
    // only one singular prime
    t0:=Cputime();
        assert #WKICM(R) eq 173;
    t_curr:=Cputime(t0);
    t_prev_best:=26;
    "Current running time: ",t_curr;
    if Abs((t_curr - t_prev_best)/t_prev_best) gt 0.1 then 
        print "The current timing is different from the previous best known one. UPDATE!"; 
    end if;
    //"the size of the output, 173 classes, has been computed using the OLD method, in 189000 seconds";

    "-------------------------------------------------------------";
    "An example in which the order has 2 singular primes";

    // Here: R has 2 singular primes
    //SetVerbose("WKICM",2);
    P<x>:=PolynomialRing(Integers());
    f:=x^10 - x^9 + 4*x^8 - 6*x^7 + 8*x^6 - 16*x^5 + 16*x^4 - 24*x^3 + 32*x^2 - 16*x + 32;
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    q:=Integers() ! Round(ConstantCoefficient(f)^(2/Degree(f)));
    R:=Order([F,q/F]);
    t0:=Cputime();
        assert #WKICM(R) eq 238;
    t_curr:=Cputime(t0);
    "Current running time: ",t_curr;
    t_prev_best:=10.4;
    if Abs((t_curr - t_prev_best)/t_prev_best) gt 0.1 then 
        print "The current timing is different from the previous best known one. UPDATE!"; 
    end if;

    "-------------------------------------------------------------";
    "An example in which the order has 3 singular primes";

    // Here: R has 3 singular primes
    //SetVerbose("WKICM",2);
    P<x>:=PolynomialRing(Integers());
    f:=x^6 - 3*x^5 + 30*x^4 - 175*x^3 + 750*x^2 - 1875*x + 15625;
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    q:=Integers() ! Round(ConstantCoefficient(f)^(2/Degree(f)));
    R:=Order([F,q/F]);
    t0:=Cputime();
        assert #WKICM(R) eq 315;
    t_curr:=Cputime(t0);
    "Current running time: ",t_curr;
    t_prev_best:=14.1;
    if Abs((t_curr - t_prev_best)/t_prev_best) gt 0.1 then 
        print "The current timing is different from the previous best known one. UPDATE!"; 
    end if;

    "-------------------------------------------------------------";
    "An example in which the order has 5 singular primes";

    // Here: R has 5 singular primes
    //SetVerbose("WKICM",2);
    P<x>:=PolynomialRing(Integers());
    f:=x^6 + 8*x^5 + 50*x^4 + 200*x^3 + 1250*x^2 + 5000*x + 15625;
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    R:=Order([F,25/F]);
    t0:=Cputime();
        assert #WKICM(R) eq 560;
    t_curr:=Cputime(t0);
    "Current running time: ",t_curr;
    t_prev_best:=16.5;
    if Abs((t_curr - t_prev_best)/t_prev_best) gt 0.1 then 
        print "The current timing is different from the previous best known one. UPDATE!"; 
    end if;

    "-------------------------------------------------------------";
    "A very big example. It should take approx 1.5h";
    P<x>:=PolynomialRing(Integers());
    f:=x^8 - 2*x^7 + 7*x^6 - 14*x^5 + 40*x^4 - 56*x^3 + 112*x^2 - 128*x + 256;
    A:=EtaleAlgebra(f);
    R:=EquationOrder(A);
    t0:=Cputime();
        assert #WKICM(R) eq 114492 ;
    t_curr:=Cputime(t0);
    "Current running time: ",t_curr;
    t_prev_best:=4462;
    if Abs((t_curr - t_prev_best)/t_prev_best) gt 0.1 then 
        print "The current timing is different from the previous best known one. UPDATE!"; 
    end if;

    "-------------------------------------------------------------";
    printf "The whole set of slow tests required %o\n",Cputime(time_start_slow_tests);
