/* vim: set syntax=magma :*/

//freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose AlgEtQTraceNorm, 3;

/*TODO:

*/

import "Ord.m" : crQZ , crZQ , Columns , hnf , MatrixAtoQ , MatrixAtoZ , MatrixQtoA , meet_zbasis , inclusion_matrix;

//------------
// Trace and Norm
//------------

intrinsic Trace(x::AlgEtQElt) -> Any
{Returns the trace of the element x of an étale algebra.}
    require HasBaseField(Algebra(x)) : "The numeber fields are not all defined over the same BaseField.";
    return &+[Trace(y) : y in Components(x)];
end intrinsic;

intrinsic Norm(x::AlgEtQElt) -> Any
{Returns the norm of the element x of an étale algebra.}
    require HasBaseField(Algebra(x)) : "The numeber fields are not all defined over the same BaseField.";
    return &*[Norm(y) : y in Components(x)];
end intrinsic;

intrinsic AbsoluteTrace(x::AlgEtQElt) -> Any
{Returns the absolute trace of the element x of an étale algebra.}
    return &+[AbsoluteTrace(y) : y in Components(x)];
end intrinsic;

intrinsic AbsoluteNorm(x::AlgEtQElt) -> Any
{Returns the absolute norm of the element x of an étale algebra.}
    return &*[AbsoluteNorm(y) : y in Components(x)];
end intrinsic;


//------------
// Trace dual ideal
//------------

intrinsic TraceDualIdeal(I::AlgEtQIdl) -> AlgEtQIdl
{Returns the trace dual ideal of an ideal in an order in an etale algebra.}
    if not assigned I`TraceDualIdeal then
        A:=Algebra(I);
        require PrimeField(A) eq BaseField(A) : "implementend only for algebras over the prime field";
        S:=Order(I);
        B:=ZBasis(I);
        Nnf:=#Components(A);
        n:=#B;
        Q:=MatrixRing(RationalField(), n)![AbsoluteTrace(B[i]*B[j]): i, j in [1..n] ];
        QQ:=Q^-1;
        //BB:=[A ! (&+[ (QQ[i,j]*B[j]): j in [1..n]]) : i in [1..n]] ; //too many coercions
        B_comp:=[Components(b) : b in B];
        BB:=< [(&+[ (QQ[i,j]*B_comp[j][k]): j in [1..n]]) : i in [1..n]] : k in [1..Nnf]>;
        BB:=[ A ! < BB[k][i] : k in [1..Nnf] > : i in [1..n] ];
        assert2 forall{ i : i,j in [1..n] | AbsoluteTrace( B[i]*BB[j] ) eq KroneckerDelta(i,j) };
        It:=Ideal(S,BB);
        It`ZBasis:=BB; //we know that BB is a ZBasis
        if assigned I`MultiplicatorRing then
        // the multiplicator ring of I is the same of its trace dual.
            It`MultiplicatorRing:=I`MultiplicatorRing;
        end if;
        ZBasisLLL(It);
        I`TraceDualIdeal:=It;
    end if;
    return I`TraceDualIdeal;
end intrinsic;

intrinsic TraceDualIdeal(O::AlgEtQOrd) -> AlgEtQIdl
{Returns the trace dual ideal of an order in an etale algebra.}
    if not assigned O`TraceDualIdeal then
        Ot:=TraceDualIdeal(OneIdeal(O));
        O`TraceDualIdeal := Ot;
    end if;
    return O`TraceDualIdeal;
end intrinsic;

/* TEST

    printf "### Testing Trace and Norm:";
    AttachSpec("~/packages_github/AlgEt/spec");
    SetVerbose("AlgEtQTraceNorm",1);

    _<x>:=PolynomialRing(Integers());
    f:=(x^8+16)*(x^8+81);
    A:=EtaleAlgebra(f);
    for i in [1..100] do
        a:=Random(A);
        b:=Random(A);
        assert Trace(a)+Trace(b) eq Trace(a+b);
        assert Norm(a)*Norm(b) eq Norm(a*b);
    end for;
    printf " all good!\n"; 

*/