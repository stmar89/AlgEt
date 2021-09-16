/* vim: set syntax=magma :*/

//freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose AlgEtTraceNorm, 3;

/*TODO:

*/

//------------
// Trace and Norm
//------------

intrinsic Trace(x::AlgEtElt) -> Any
{Returns the trace of the element x of an étale algebra.}
    require HasBaseField(Algebra(x)) : "The numeber fields are not all defined over the same BaseField.";
    return &+[Trace(y) : y in Components(x)];
end intrinsic;

intrinsic Norm(x::AlgEtElt) -> Any
{Returns the norm of the element x of an étale algebra.}
    require HasBaseField(Algebra(x)) : "The numeber fields are not all defined over the same BaseField.";
    return &*[Norm(y) : y in Components(x)];
end intrinsic;

intrinsic AbsoluteTrace(x::AlgEtElt) -> Any
{Returns the absolute trace of the element x of an étale algebra.}
    return &+[AbsoluteTrace(y) : y in Components(x)];
end intrinsic;

intrinsic AbsoluteNorm(x::AlgEtElt) -> Any
{Returns the absolute norm of the element x of an étale algebra.}
    return &*[AbsoluteNorm(y) : y in Components(x)];
end intrinsic;


/* TEST

    Attach("~/packages_github/AlgEt/AlgEt.m");
    Attach("~/packages_github/AlgEt/AlgEtElt.m");
    Attach("~/packages_github/AlgEt/AlgEtOrd.m");
    Attach("~/packages_github/AlgEt/AlgEtTraceNorm.m");
    SetVerbose("AlgEtTraceNorm",2);

    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    for i in [1..100] do
        a:=Random(A);
        b:=Random(A);
        assert Trace(a)+Trace(b) eq Trace(a+b);
        assert Norm(a)*Norm(b) eq Norm(a*b);
    end for;

    K:=NumberField(x^2-5);
    _<y>:=PolynomialRing(K);
    E1:=NumberField(y^2-49*7*K.1);
    E2:=NumberField(y^5-25*7*K.1);
    A:=EtaleAlgebra([E1,E2]); 
    assert HasBaseField(A);
    O:=MaximalOrder(A);
    for i in [1..100] do
        a:=Random(O);
        b:=Random(O);
        assert AbsoluteTrace(a)+AbsoluteTrace(b) eq AbsoluteTrace(a+b);
        assert Trace(a)+Trace(b) eq Trace(a+b);
        assert AbsoluteNorm(a)*AbsoluteNorm(b) eq AbsoluteNorm(a*b);
        assert Norm(a)*Norm(b) eq Norm(a*b);
    end for;
    
    A:=EtaleAlgebra([K,K]);
    for i in [1..100] do
        a:=Random(A);
        b:=Random(A);
        assert Trace(a)+Trace(b) eq Trace(a+b);
        assert Norm(a)*Norm(b) eq Norm(a*b);
    end for;


    seq:=[NumberField(p),NumberField(x^2-5)];
    A:=EtaleAlgebra(seq);
    time O2:=Order(AbsoluteBasis(A));
    for O in [O2] do
        for i in [1..100] do
            assert Random(O) in O;
        end for;
    end for;

    K:=NumberField(x^2-25*5);
    _<y>:=PolynomialRing(K);
    p:=y^2-7;
    E:=NumberField(p); //relative ext
    A:=EtaleAlgebra([E,E]);
    O:=MaximalOrder(A);
    for i in [1..100] do
        a:=Random(O);
        b:=Random(O);
        assert AbsoluteTrace(a)+AbsoluteTrace(b) eq AbsoluteTrace(a+b);
        assert Trace(a)+Trace(b) eq Trace(a+b);
        assert AbsoluteNorm(a)*AbsoluteNorm(b) eq AbsoluteNorm(a*b);
        assert Norm(a)*Norm(b) eq Norm(a*b);
    end for;

*/
