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

declare verbose GraphOverOrders, 1;

declare attributes AlgEtQOrd : GraphOverOrders;

//------------
// GraphOverOrders : the graph of (minimal) inclusions of overorders of R
//------------

intrinsic GraphOverOrders(R:AlgEtQOrd) -> GrphDir
{Given an order R returns the graph G of minimal inclusions of the overorders of R. More precisely, the vertices of G are integers between 1 and the number of OverOrders(R), and there is an edge [i,j] if and only if OverOrder(R)[j] is a minimal overorder of OverOrders(R)[i].}
    if not assigned R`GraphOverOrders then
        if IsMaximal(R) then
            G:=Digraph< 1 | >;
        else
            pp:=SingularPrimes(R);
            oo_at_Ps:=[ OverOrdersAtPrime(R,P) : P in pp ]; // 8 secs
            assert forall{ x : x in oo_at_Ps | x[1] eq R};
            edges:=[];
            for iP->P in pp do
                oo_P:=oo_at_Ps[iP];
                edges_P:=[ ]; // minimal P-inclusions
                for i->S in oo_P do
                    if assigned S`MinimalOverOrders then
                        PS:=S!!P;
                        min:=&join[ min[2] : min in S`MinimalOverOrders | PS subset min[1] ];
                        for T in min do
                            j:=Index(oo_P,T);
                            assert j ne 0;
                            Append(~edges_P,[i,j]);
                        end for;
                    end if;
                end for;
                // test
                assert2 forall{ e : e in edges_P | oo_P[e[1]] subset oo_P[e[2]]}; // inclusions
                // the next test is very time consuming
                assert2 forall{ e : e in edges_P | [ S : S in oo_P | S ne oo_P[e[1]] and oo_P[e[1]] subset S 
                                and S subset oo_P[e[2]] ] eq [oo_P[e[2]]]}; // minimal inclusions
                Append(~edges,edges_P);
            end for;
            // the following is taken from OverOrders()
            output:=OverOrders(R);
            if #pp eq 1 then
                min_inclusions:=edges[1];
                //test
                assert forall{i : i in [1..#output] | output[i] eq oo_at_Ps[1][i]};
            else
                cc:=CartesianProduct([[1..#x] : x in oo_at_Ps]); // sorted as output (check for O!)
                // given c1,c2 in cc we have that c2 is a minimaloverorder of c1 if and only if 
                // - c1[i] eq c2[i] for all i's but one, say i0, and
                // - c2[i0] is a minimaloverorder of c1[i0],
                // which is equivalent to have [c1[i0],c2[i0]] in edges[i0].
                min_inclusions:=[];
                cc:=[ c : c in cc ];
                for i,j in [1..#output] do
                    c1:=cc[i];
                    c2:=cc[j];
                    different_entries:=[ k : k in [1..#c1] | c1[k] ne c2[k] ];
                    if #different_entries eq 1 then
                        // c1 and c2 differ only at the prime P=pp[i0] of R.
                        i0:=different_entries[1];
                        edges_P:=edges[i0];
                        incl_atP:=[c1[i0],c2[i0]];
                        if incl_atP in edges[i0] then
                            // we found a minimal inclusion: output[i] subset output[j]]
                            Append(~min_inclusions,[i,j]);
                        end if;
                    end if;
                end for;
                // test
                assert2 forall{ incl : incl in min_inclusions | output[incl[1]] subset output[incl[2]]}; // inclusion
                // the next test is very time consuming
                assert2 forall{ incl : incl in min_inclusions | [ S : S in output | 
                                S ne output[incl[1]] and output[incl[1]] subset S 
                                    and S subset output[incl[2]]] eq [output[incl[2]]]}; // minimal inclusions
            end if;
            G:=Digraph< #output | min_inclusions >;
        end if;
        R`GraphOverOrders:=G;
    end if;
    return R`GraphOverOrders;
end intrinsic;

/* TESTS
    
    printf "### Testing GraphOverOrders:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    P<x>:=PolynomialRing(Integers());
    fs:=[ 
          x^8 + 16, //1 sing prime
          x^4-10000*x^3-10000*x^2-10000*x-10000, //2 sing primes
          x^4-30^3*x^3-30^3*x^2-30^3*x-30^3 //3 sing primes
        ];
    for f in fs do
        printf ".";
        A:=EtaleAlgebra(f);
        _:=GraphOverOrders(MaximalOrder(A));
        R:=EquationOrder(A); 
        _:=OverOrders(R);
        SetAssertions(2);
        _:=GraphOverOrders(R);
        SetAssertions(1);
    end for;
    printf " all good!\n";

*/
