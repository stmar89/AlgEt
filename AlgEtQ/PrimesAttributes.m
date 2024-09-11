/* vim: set syntax=magma : */

/////////////////////////////////////////////////////
// Stefano Marseglia, stefano.marseglia89@gmail.com
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
// Copyright 2024, S. Marseglia
/////////////////////////////////////////////////////

declare attributes AlgEtQIdl : InertiaDegree,
                               RamificationIndex;

declare verbose nice_unif,3;

intrinsic Valuation(x::AlgEtQElt,P::AlgEtQIdl)->RngIntElt
{Given an element x and a maximal ideal P of the maximal order, it returns the valuation of x at P.}
    OO:=Order(P);
    require IsMaximal(OO) and IsPrime(P): "The input needs to be a maximal ideal of the maximal order.";
    _,comp:=IsProductOfIdeals(P);
    ind:=[ i : i in [1..#comp] | not One(NumberField(Order(comp[i]))) in comp[i] ];
    assert #ind eq 1;
    ind:=ind[1];
    return Valuation(Components(x)[ind],comp[ind]);
end intrinsic;

intrinsic Valuation(I::AlgEtQIdl,P::AlgEtQIdl)->RngIntElt
{Given a fractional ideal I and a maximal ideal P, both of the maximal order, it returns the valuation of I at P.}
    OO:=Order(P);
    require IsMaximal(OO) and IsPrime(P): "The input needs to be a maximal ideal of the maximal order.";
    _,compI:=IsProductOfIdeals(I);
    _,compP:=IsProductOfIdeals(P);
    ind:=[ i : i in [1..#compP] | not One(NumberField(Order(compP[i]))) in compP[i] ];
    assert #ind eq 1;
    ind:=ind[1];
    return Valuation(compI[ind],compP[ind]);
end intrinsic;

intrinsic InertiaDegree(P::AlgEtQIdl)->RngIntElt
{Given a maximal ideal P of the maximal order O above the rational prime p, it returns the inertia degree of P, that is, the index of the finite field extension GF(p)->O/P.}
    if not assigned P`InertiaDegree then
        OO:=Order(P);
        require IsMaximal(OO) and IsPrime(P): "The input needs to be a maximal ideal of the maximal order.";
        q:=Index(OO,P);
        t,p,f:=IsPrimePower(q);
        assert t;
        P`InertiaDegree:= f;
    end if;
    return P`InertiaDegree;
end intrinsic;

intrinsic RamificationIndex(P::AlgEtQIdl)->RngIntElt
{Given a maximal ideal P of the maximal order O, it returns the reamification index of P.}
    if not assigned P`RamificationIndex then
        OO:=Order(P);
        require IsMaximal(OO) and IsPrime(P): "The input needs to be a maximal ideal of the maximal order.";
        q:=Index(OO,P);
        t,p:=IsPrimePower(q);
        assert t;
        P`RamificationIndex:=Valuation(Algebra(P)!p,P);
    end if;
    return P`RamificationIndex;
end intrinsic;

intrinsic Uniformizers(PPs::SeqEnum[AlgEtQIdl])->SeqEnum
{Given a sequence of maximal ideals P of the maximal order, it returns a sequence of elements t_P such that t_P is a uniformizer of P and a unit at every other prime in the sequence.}
    OO:=Order(PPs[1]);
    require IsMaximal(OO) and forall{P : P in PPs | IsPrime(P)} : "The input needs to be a sequenc of maximal ideals of the maximal order.";
    vprintf nice_unif,2 : "nice_uniformizers: %o primes\n",#PPs;
    one:=One(Algebra(PPs[1]));
    nice_unifs:=[];
    PPs2:=[ P^2 : P in PPs];
    for iP->P in PPs do
        vprintf nice_unif,2 : "nice_uniformizers: %o-th prime...",iP;
        P2:=PPs2[iP];
        Q,q:=Quotient(P,P2);
        repeat
            t:=Random(Q);
        until t ne Zero(Q);
        vprintf nice_unif,2 : "got non-zero elt...";
        tP:=t@@q;
        while IsZeroDivisor(tP) do
            tP +:=Random(P2);
        end while;
        vprintf nice_unif,2 : "got non-zerodivisor...";
        assert2 tP in P and not tP in P2;
        elts:=[ i eq iP select tP else one : i in [1..#PPs] ];
        tP:=CRT(PPs2,elts);
        vprintf nice_unif,2 : "asserts";
        assert2 tP in P and not tP in P2;
        vprintf nice_unif,2 : ".";
        assert2 tP in MaximalOrder(Algebra(tP));
        vprintf nice_unif,2 : ".";
        assert2 forall{ Q : Q in [ Q : Q in PPs | Q ne P ] | not tP in Q };
        vprintf nice_unif,2 : ".";
        assert3 Seqset(PrimesAbove(tP*MaximalOrder(Algebra(tP)))) meet Seqset(PPs) eq {P}; //for some reason this is very time consuming
        vprintf nice_unif,2 : ".";
        assert not IsZeroDivisor(tP);
        vprintf nice_unif,2 : "all passed\n";
        Append(~nice_unifs,tP);
    end for;
    return nice_unifs;
end intrinsic;


