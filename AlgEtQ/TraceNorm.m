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

freeze;

declare verbose AlgEtQTraceNorm, 3;

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
{Returns the absolute trace of the element x of an étale algebra. Since the étale algebra is over the rationals this is the same as Trace.}
    return &+[AbsoluteTrace(y) : y in Components(x)];
end intrinsic;

intrinsic AbsoluteNorm(x::AlgEtQElt) -> Any
{Returns the absolute norm of the element x of an étale algebra. Since the étale algebra is over the rationals this is the same as Norm.}
    return &*[AbsoluteNorm(y) : y in Components(x)];
end intrinsic;


//------------
// Trace dual ideal
//------------

intrinsic TraceDualIdeal(I::AlgEtQIdl) -> AlgEtQIdl
{Returns the trace dual ideal of the ideal I, that is, the set of elements x of the algebra such that Trace(x*I) is integer-valued.}
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
{Returns the trace dual ideal of an order in an etale algebra, that is, the set of elements x of the algebra such that Trace(x*O) is integer-valued.}
    if not assigned O`TraceDualIdeal then
        Ot:=TraceDualIdeal(OneIdeal(O));
        O`TraceDualIdeal := Ot;
    end if;
    return O`TraceDualIdeal;
end intrinsic;

/* TESTS

    printf "### Testing Trace and Norm:";
    //AttachSpec("~/packages_github/AlgEt/spec");
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
