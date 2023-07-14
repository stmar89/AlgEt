/* vim: set syntax=magma :*/

    //SetAssertions(2);

    //AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/dev/rec_wk_icm.m");
    //SetVerbose("WKICM",2);
    P<x>:=PolynomialRing(Integers());
    f:=x^8+16;
    A:=EtaleAlgebra(f);
    R:=EquationOrder(A);
    // only one singular prime
    t0:=Cputime();
        assert #WKICM(R) eq 173;
    t_curr:=Cputime(t0);
    t_prev_best:=22.5;
    "Current running time: ",t_curr;
    if Abs((t_curr - t_prev_best)/t_prev_best) gt 0.1 then 
        print "The current timing is different from the previous best known one. UPDATE!"; 
    end if;
    //"the size of the output, 173 classes, has been computed using the OLD method, in 189000 seconds";

    // Here: R has 2 singular primes
    //AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/dev/rec_wk_icm.m");
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

    // Here: R has 3 singular primes
    //AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/dev/rec_wk_icm.m");
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
    t_prev_best:=28.4;
    if Abs((t_curr - t_prev_best)/t_prev_best) gt 0.1 then 
        print "The current timing is different from the previous best known one. UPDATE!"; 
    end if;

    // Here: R has 5 singular primes
    //AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/dev/rec_wk_icm.m");
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
    t_prev_best:=27;
    if Abs((t_curr - t_prev_best)/t_prev_best) gt 0.1 then 
        print "The current timing is different from the previous best known one. UPDATE!"; 
    end if;

/* Running on screen big_test_old_rec using rec_wk_icm.m
   Update timings when that one is done
    // VERY big it should not be in this test suit
    // I have improved a lot KnownOrders. I expect that populating OverORders should be much faster now
    //AttachSpec("~/packages_github/AlgEt/spec");
    Attach("~/packages_github/AlgEt/dev/rec_wk_icm.m");
    //SetVerbose("WKICM",2);
    P<x>:=PolynomialRing(Integers());
    f:=x^8 - 2*x^7 + 7*x^6 - 14*x^5 + 40*x^4 - 56*x^3 + 112*x^2 - 128*x + 256;
    A:=EtaleAlgebra(f);
    R:=EquationOrder(A);
    time assert #WKICM(R) eq 114492; // 114492 in ~519000 second ~ 144h with the OLD method
                                     // Now we have 1.63 h
*/

