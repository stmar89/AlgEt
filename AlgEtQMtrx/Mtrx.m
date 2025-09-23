/* vim: set syntax=magma :*/

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

freeze;

/*
    Write the étale algebra as K=K1x...xKn with each Ki a number field.
    Then each sxr matrix M with extries in K is a 'vector' of matrices (M1,...,Mn) where 
    each Mi is a sxr matrix with entries in Ki. 
    Operations are then performed component-wise.
*/

declare verbose AlgEtQMtrx, 3;

declare type AlgEtQMtrx;

declare attributes AlgEtQMtrx : 
                               Universe, // AlgEtQ
                               NumberOfRows,
                               NumberOfColumns,
                               Components,
                               Entries, // SeqEnum[SeqEnum[AlgEtQlt]]
                               Rank,
                               Trace,
                               Determinant;
                               // MinimalPolynomial does not work well if there are zerodivisors.
                               // CharacteristicPolynomial require Polynomials over AlgEtQ, not implemented.

//------------
// Printing 
//------------

intrinsic Print(x::AlgEtQMtrx)
{Print the element.}
    /*
    if not assigned x`Entries then
        assert assigned x`Components;
        print Components(x);
    end if;
    */
    Q:=Entries(x);
    Q:=[ [ &cat(Split(Sprint(e)," ")) : e in seq ] : seq in Q ];
    max:=Max([ #s : s in Flat(Q) ]);
    out:="";
    for row in Q do
        out cat:="[";
        for e in row do
             out cat:=&cat([""] cat [ " " : i in [1..max-#e]]) cat e cat " ";
        end for;
        Prune(~out); // remove last space
        out cat:="]\n";
    end for;
    Prune(~out); //remove he last \n
    printf "%o",out;
end intrinsic;

//------------
// Access attributes
//------------

intrinsic Universe(x::AlgEtQMtrx) -> AlgEtQ
{Returns the étale algebra to which the entries of the matrix belongs to.}
  return x`Universe;
end intrinsic;

intrinsic NumberOfRows(x::AlgEtQMtrx) -> RngIntElt
{Returns the number of rows of the matrix.}
  return x`NumberOfRows;
end intrinsic;

intrinsic NumberOfColumns(x::AlgEtQMtrx) -> RngIntElt
{Returns the number of colums of the matrix.}
  return x`NumberOfColumns;
end intrinsic;

intrinsic Components(x::AlgEtQMtrx) -> Tup
{Returns the tuple of components of x. If x is defined over K=K1x...xKn with each Ki a number field, the ith compoenent is a matrix over Ki.}
    if not assigned x`Components then
        assert assigned x`Entries and assigned x`Universe;
        entries:=Entries(x);
        n:=#Components(Universe(x));
        Q:=[ Components(e) : e in Flat(entries)];
        Q:=< [ e[i] : e in Q ] : i in [1..n] >;
        Q:=< Matrix(#entries,#entries[1],Q[i]) : i in [1..n] >;
        x`Components:=Q;
    end if;
    return x`Components;
end intrinsic;

intrinsic Entries(x::AlgEtQMtrx) -> SeqEnum[SeqEnum[AlgEtElt]]
{Returns a sequence of sequences containing the entries of the matrix.}
    if not assigned x`Entries then
        assert assigned x`Components and assigned x`Universe;
        comp:=Components(x);
        A:=Universe(x);
        x`Entries:=[ [A!<c[i,j]:c in comp> : j in [1..NumberOfColumns(comp[1])]] : i in [1..NumberOfRows(comp[1])] ];
    end if;
    return x`Entries;
end intrinsic;

//------------
// Creation
//------------

intrinsic Matrix(entries::SeqEnum[SeqEnum[AlgEtQElt]]) -> AlgEtQMtrx
{Given a sequence of s sequences of r elements of an étale algebra, it returns the corresponding sxr matrix.}
    require #entries gt 0 : "Empty sequence of entries.";
    require forall{ s : s in entries | #s eq #entries[1] }: "Each sequence in entries must have the same length.";
    A:=Algebra(entries[1,1]);
    x:=New(AlgEtQMtrx);
    x`Universe:=A;
    x`NumberOfRows:=#entries;;
    x`NumberOfColumns:=#entries[1];
    x`Entries:=entries;
    return x;
end intrinsic;

intrinsic Matrix(s::RngIntElt, r::RngIntElt, Q::SeqEnum[AlgEtQElt]) -> AlgEtQMtrx
{Given a sequence of sxr elements of an étale algebra, it returns the corresponding sxr matrix.}
    require s*r eq #Q: "Not enough entries.";
    Q:=Partition(Q,r);
    assert #Q eq s;
    return Matrix(Q);
end intrinsic;

CreateMtrxFromComponents:=function(A,comp)
// as internal, since giving a Tup an input signature will probably create problems
// A is the algebra and comp the Componens
    s:=NumberOfRows(comp[1]);
    r:=NumberOfColumns(comp[1]);
    assert forall{ c : c in comp | NumberOfRows(c) eq s and NumberOfColumns(c) eq r };
    assert Type(A) eq AlgEtQ;
    x:=New(AlgEtQMtrx);
    x`Universe:=A;
    x`NumberOfRows:=s;
    x`NumberOfColumns:=r;
    x`Components:=comp;
    return x;
end function;

intrinsic DiagonalMatrix(Q::SeqEnum[AlgEtQElt]) -> AlgEtQMtrx
{Given a sequence of elements of an étale algebra, it returns the corresponding diagonal matrix.}
    A:=Algebra(Q[1]);
    n:=#Components(A);
    diags:=< [ Components(d)[i] : d in Q ] : i in [1..n] >;
    comp:=< DiagonalMatrix(diags[i]) : i in [1..n] >;
    return CreateMtrxFromComponents(A,comp);
end intrinsic;

//------------
// Access, Rows, Columns, Transpose
//------------

/*
TODO How to???
intrinsic '[,]'(x::AlgEtQMtrx,i::RngIntElt,j::RngIntElt) -> AlgEtQElt
{Returns x[i,j].}
    return Entries(x)[i,j];
end intrinsic;
*/

intrinsic Rows(x::AlgEtQMtrx) -> SeqEnum[SeqEnum[AlgEtQElt]]
{Returns the sequence of rows.}
    return Entries(x);
end intrinsic;

intrinsic Columns(x::AlgEtQMtrx) -> SeqEnum[SeqEnum[AlgEtQElt]]
{Returns the sequence of columns.}
    en:=Entries(x);
    en:=[ [ en[i,j] : i in [1..NumberOfRows(x)] ] : j in [1..NumberOfColumns(x)] ]; 
    return en;
end intrinsic;

intrinsic Transpose(x::AlgEtQMtrx) -> AlgEtQMtrx
{Returns the transpose matrix.}
    return Matrix(Columns(x));
end intrinsic;

intrinsic IsSquareMatrix(x::AlgEtQMtrx) -> BoolElt
{Returns whether the matrix is a square matrix.}
  return NumberOfRows(x) eq NumberOfColumns(x);
end intrinsic;

//------------
// Basic and Random elements
//------------

intrinsic ScalarMatrix(n::RngIntElt,s::AlgEtQElt) -> AlgEtQMtrx
{Returns the scalar matrix s*I where I is the nxn identity matrix.}
    Q:=[ s : i in [1..n] ];
    return DiagonalMatrix(Q);
end intrinsic;

intrinsic IdentityMatrix(A::AlgEtQ, n::RngIntElt) -> AlgEtQMtrx
{Returns the nxn identity matrix over the étale algebra A.}   
    u:=One(A);
    return ScalarMatrix(n,u);
end intrinsic;

intrinsic ZeroMatrix(A::AlgEtQ, n::RngIntElt) -> AlgEtQMtrx
{Returns the nxn identity matrix over the étale algebra A.}   
    u:=Zero(A);
    return ScalarMatrix(n,u);
end intrinsic;

intrinsic RandomMatrix(A::AlgEtQ, s::RngIntElt, r::RngIntElt : bd:=3) -> AlgEtQMtrx
{Returns a random sxr identity matrix over the étale algebra A, with coefficients bounded by the VarArg bd.}   
    Q:=[ [ Random(A,bd) : j in [1..r] ] : i in [1..s] ];
    return Matrix(Q);
end intrinsic;

intrinsic RandomMatrix(A::AlgEtQ, n::RngIntElt : bd:=3) -> AlgEtQMtrx
{Returns a random square nxn identity matrix over the étale algebra A, with coefficients bounded by the VarArg bd.}   
    return RandomMatrix(A,n,n : bd:=bd);
end intrinsic;

//------------
// IsInvertble, Rank, Trace and Determinant
//------------
intrinsic IsInvertible(x::AlgEtQMtrx) -> BoolElt,AlgEtQMtrx
{Returns whether the matrix is invertible and, if so, it returns also the inverse.}
    if not IsSquareMatrix(x) then
        return false,_;
    end if;
    A:=Universe(x);
    x:=Components(x);
    x_inv:=<>;
    for c in x do
        test,c_inv:=IsInvertible(c);
        if not test then
            return false,_;
        else
            Append(~x_inv,c_inv);
        end if;
    end for;
    // if we exit the loop then all components are invertible
    return true,CreateMtrxFromComponents(A,x_inv);
end intrinsic;

intrinsic Rank(x::AlgEtQMtrx) -> RngIntElt
{The rank of x, which is defined as the minimum of the rank of his components.}
    m:=Min([Rank(c) : c in Components(x)]);
    return m;
end intrinsic;

intrinsic Trace(x::AlgEtQMtrx) -> AlgEtQElt
{Determinant.}
    A:=Universe(x);
    return A!<Trace(c) : c in Components(x)>;
end intrinsic;

intrinsic Determinant(x::AlgEtQMtrx) -> AlgEtQElt
{Determinant.}
    A:=Universe(x);
    return A!<Determinant(c) : c in Components(x)>;
end intrinsic;

//------------
// Equality and Operations
//------------

// eq 
intrinsic 'eq'(x1::AlgEtQMtrx,x2::AlgEtQMtrx) -> BoolElt
{Is x1=x2 ?}
    A:=Universe(x1);
    require A cmpeq Universe(x2): "The matrices are not defined over the same algebra.";
    if assigned x1`Components and assigned x2`Components then // matrix operations are only on components.
                                                              // it is more likely that entries are not assigned.
        return Components(x1) eq Components(x2);
    else
        return Entries(x1) eq Entries(x2);
    end if;
end intrinsic;

// // +
// // using the components, as below, is faster!
// intrinsic '+'(x1::AlgEtQMtrx,x2::AlgEtQMtrx) -> AlgEtQMtrx
// {x1+x2.}
//     A:=Universe(x1);
//     require A cmpeq Universe(x2): "The matrices are not defined over the same algebra.";
//     s:=NumberOfRows(x1);
//     r:=NumberOfColumns(x1);
//     sr:=s*r;
//     require NumberOfRows(x2) eq s and NumberOfColumns(x2) eq r : "The matrices cannot be added.";
//     x1:=Flat(Entries(x1));
//     x2:=Flat(Entries(x2));
//     x3:=[ x1[k]+x2[k] : k in [1..sr] ]; 
//     return Matrix(s,r,x3);
// end intrinsic;
// 
// // -
// intrinsic '-'(x1::AlgEtQMtrx) -> AlgEtQMtrx
// {-x1.}
//     A:=Universe(x1);
//     s:=NumberOfRows(x1);
//     r:=NumberOfColumns(x1);
//     sr:=s*r;
//     x1:=Flat(Entries(x1));
//     x3:=[ -x1[k] : k in [1..sr] ]; 
//     return Matrix(s,r,x3);
// end intrinsic;
// 
// intrinsic '-'(x1::AlgEtQMtrx,x2::AlgEtQMtrx) -> AlgEtQMtrx
// {x1-x2.}
//     A:=Universe(x1);
//     require A cmpeq Universe(x2): "The matrices are not defined over the same algebra.";
//     s:=NumberOfRows(x1);
//     r:=NumberOfColumns(x1);
//     sr:=s*r;
//     require NumberOfRows(x2) eq s and NumberOfColumns(x2) eq r : "The matrices cannot be added.";
//     x1:=Flat(Entries(x1));
//     x2:=Flat(Entries(x2));
//     x3:=[ x1[k]-x2[k] : k in [1..sr] ]; 
//     return Matrix(s,r,x3);
// end intrinsic;
// 
// // *
// intrinsic '*'(x1::AlgEtQMtrx,x2::AlgEtQMtrx) -> AlgEtQMtrx
// {x1*x2.}
//     A:=Universe(x1);
//     require A cmpeq Universe(x2): "The elements must belong to the same algebra.";
//     require NumberOfColumns(x1) eq NumberOfRows(x2) : "The matrices have the wrong sizes.";
//     s:=NumberOfRows(x1);
//     r:=NumberOfColumns(x2);
//     rows1:=Rows(x1);
//     cols2:=Columns(x2);
//     x3:=[ [ DotProduct(rows1[i],cols2[j]) : j in [1..r] ] : i in [1..s] ];
//     x3:=Matrix(x3);
//     assert2 NumberOfRows(x3) eq s and NumberOfColumns(x3) eq r;
//     return x3;
// end intrinsic;

// +
intrinsic '+'(x1::AlgEtQMtrx,x2::AlgEtQMtrx) -> AlgEtQMtrx
{x1+x2.}
    A:=Universe(x1);
    require A cmpeq Universe(x2): "The matrices are not defined over the same algebra.";
    comp1:=Components(x1);
    comp2:=Components(x2);
    return CreateMtrxFromComponents(A,<comp1[i]+comp2[i] : i in [1..#comp1]> );
end intrinsic;

// -
intrinsic '-'(x1::AlgEtQMtrx) -> AlgEtQMtrx
{-x1.}
    A:=Universe(x1);
    return CreateMtrxFromComponents(A,<-c : c in Components(x1)>);
end intrinsic;

intrinsic '-'(x1::AlgEtQMtrx,x2::AlgEtQMtrx) -> AlgEtQMtrx
{x1-x2.}
    A:=Universe(x1);
    require A cmpeq Universe(x2): "The matrices are not defined over the same algebra.";
    comp1:=Components(x1);
    comp2:=Components(x2);
    return CreateMtrxFromComponents(A,<comp1[i]-comp2[i] : i in [1..#comp1]> );
end intrinsic;

// *
intrinsic '*'(x1::AlgEtQMtrx,x2::AlgEtQMtrx) -> AlgEtQMtrx
{x1*x2.}
    A:=Universe(x1);
    require A cmpeq Universe(x2): "The matrices are not defined over the same algebra.";
    comp1:=Components(x1);
    comp2:=Components(x2);
    return CreateMtrxFromComponents(A,<comp1[i]*comp2[i] : i in [1..#comp1]> );
end intrinsic;

intrinsic '*'(x1::.,x2::AlgEtQMtrx) -> AlgEtQMtrx
{x1*x2.}
    bool,x1:=IsCoercible(Universe(x2),x1);
    require bool : "x1 not coercible";
    comp1:=Components(x1);
    comp2:=Components(x2);
    return CreateMtrxFromComponents(Universe(x2),<comp1[i]*comp2[i] : i in [1..#comp1]> );
end intrinsic;

intrinsic '*'(x1::AlgEtQMtrx,x2::.) -> AlgEtQMtrx
{x1*x2.}
    bool,x2:=IsCoercible(Universe(x1),x2);
    require bool : "x2 not coercible";
    comp1:=Components(x1);
    comp2:=Components(x2);
    return CreateMtrxFromComponents(Universe(x2),<comp1[i]*comp2[i] : i in [1..#comp1]> );
end intrinsic;

intrinsic '*'(x1::RngIntElt,x2::AlgEtQMtrx) -> AlgEtQMtrx
{x1*x2.}
    return CreateMtrxFromComponents(Universe(x2),<x1*c : c in Components(x2)>);
end intrinsic;

intrinsic '*'(x1::FldRatElt,x2::AlgEtQMtrx) -> AlgEtQMtrx
{x1*x2.}
    return CreateMtrxFromComponents(Universe(x2),<x1*c : c in Components(x2)>);
end intrinsic;

intrinsic '*'(x1::AlgEtQMtrx,x2::RngIntElt) -> AlgEtQMtrx
{x1*x2.}
    return CreateMtrxFromComponents(Universe(x1),<c*x2 : c in Components(x1)>);
end intrinsic;

intrinsic '*'(x1::AlgEtQMtrx,x2::FldRatElt) -> AlgEtQMtrx
{x1*x2.}
    return CreateMtrxFromComponents(Universe(x1),<c*x2 : c in Components(x1)>);
end intrinsic;

intrinsic '^'(x::AlgEtQMtrx,n::RngIntElt) -> AlgEtQMtrx
{x^n.}
    A:=Universe(x);
    require IsSquareMatrix(x) : "The matrix is not a square.";
    if n eq 0 then
        return IdentityMatrix(A,NumberOfRows(x));
    elif n gt 0 then
        if n eq 1 then
            return x;
        end if;
        comp:=< Components(x)[i]^n : i in [1..#Components(A)] >;
        y:=CreateMtrxFromComponents(A,comp);
        return y;
    elif n lt 0 then
        is_inv,x:=IsInvertible(x);
        require is_inv : "The element is not invertible.";
        n:=-n;
        if n eq 1 then
            return x;
        end if;
        comp:=< c^n : c in Components(x) >;
        y:=CreateMtrxFromComponents(A,comp);
        return y;
    end if;
end intrinsic;

//------------
// Sums and Products of sequences
//------------

intrinsic '&+'(seq::SeqEnum[AlgEtQMtrx]) -> AlgEtQMtrx
{Given a sequence of AlgEtQMtrx returns the sum of the matrices.}
    require #seq ne 0 : "The sequence should be non-empty.";
    if #seq eq 1 then
        out:=seq[1];
    else
        comps:=<>;
        N:=#Components(seq[1]);
        comps:=< &+[ Components(s)[i] : s in seq ] : i in [1..N] >;
        out:=CreateMtrxFromComponents(Universe(seq[1]),comps);
    end if;
    return out;
end intrinsic;

intrinsic '&*'(seq::SeqEnum[AlgEtQMtrx]) -> AlgEtQMtrx
{Given a sequence of AlgEtQMtrx returns the product of the matrices.}
    require #seq ne 0 : "The sequence should be non-empty.";
    if #seq eq 1 then
        out:=seq[1];
    else
        comps:=<>;
        N:=#Components(seq[1]);
        comps:=< &*[ Components(s)[i] : s in seq ] : i in [1..N] >;
        out:=CreateMtrxFromComponents(Universe(seq[1]),comps);
    end if;
    return out;
end intrinsic;

//------------
// InsertBlock
//------------

intrinsic InsertBlock(x::AlgEtQMtrx,y::AlgEtQMtrx,i::RngIntElt,j::RngIntElt) -> AlgEtQMtrx
{Return the matrix obtained by inserting matrix Y in X at position [i,j].}
    A:=Universe(x);
    x:=Components(x);
    y:=Components(y);
    n:=#x;
    z:=<InsertBlock(x[k],y[k],i,j) : k in [1..n]>;
    return CreateMtrxFromComponents(A,z);
end intrinsic;

//------------
// Solution
//------------

intrinsic Solution(x::AlgEtQMtrx,y::AlgEtQMtrx)->AlgEtQMtrx
{Returns a solution V to the system V.x=y.}
    A:=Universe(x);
    N:=#Components(x);
    x:=Components(x);
    y:=Components(y);
    comps:=< Solution(x[i],y[i]) : i in [1..N] >;
    return CreateMtrxFromComponents(A,comps);
end intrinsic;

/*
intrinsic ???(x::AlgEtQMtrx) -> ???
{????.}
  return ????;
end intrinsic;
*/

/* TESTS

    printf "### Testing Matrices:";
    AttachSpec("~/AlgEt/spec");
    AttachSpec("~/AlgEt/specMtrx");
    SetVerbose("AlgEtQMtrx",2);

    _<x>:=PolynomialRing(Integers());
    f:=x^2+1;
    K<i>:=NumberField(f);
    g:=x^2+5;
    L<z>:=NumberField(g);
    A:=EtaleAlgebra([K,L]);
    S:=Matrix([ [A!<i,z>,Zero(A)],[Zero(A),A!<2*i,3*z>]]);
    M:=Matrix([ [A!<i,z>,Zero(A)],[Zero(A),A!<2*i,3*z>],[Zero(A),A!<1+i,1-z>]]);
    M;
    N:=Matrix(NumberOfRows(M),NumberOfColumns(M),Flat(Entries(M)));
    assert M eq N;
    assert Columns(M) eq Rows(Transpose(M));
    assert M + M eq M + N;
    OO:=ZeroMatrix(A,NumberOfRows(S));
    II:=IdentityMatrix(A,NumberOfRows(S));
    assert (S+II)-(S+ScalarMatrix(NumberOfColumns(S),One(A))) eq OO;
    assert S*OO eq OO*S and S*OO eq OO;
    assert S*II eq II*S and S*II eq S;
    _:=M*Transpose(M);

    // R:=RandomMatrix(A,100); time RR:=R+R; time _:=-R; time _:=RR-R; time _:=RR*R;
    // R:=RandomMatrix(A,200); time RR:=R+R; time _:=-R; time _:=RR-R; time _:=RR*R;
    // R:=RandomMatrix(A,300); time RR:=R+R; time _:=-R; time _:=RR-R; time _:=RR*R;
   
    R:=RandomMatrix(A,30,20);
    assert 2*R eq R*2;
    assert (-1/5)*R eq R*(1/-5);

    II:=IdentityMatrix(A,10);
    R:=RandomMatrix(A,10);
    test,Ri:=IsInvertible(R);
    assert R*Ri eq II and Ri*R eq II;
    assert R^-7 eq Ri^7;
    assert Determinant(R) eq 1/Determinant(Ri);
    assert Determinant(II) eq One(A);
    assert Determinant(OO) eq Zero(A);
    assert Trace(II) eq A!10;

    S:=Matrix([ [A!<i,z>,Zero(A)],[Zero(A),A!<2*i,0>]]);
    Rank(S);
    Determinant(S);

    assert &+[ i*R : i in [1..4] ] eq &+[1..4]*R;
    assert &*[ i*R : i in [1..4] ] eq &*[1..4]*R^4;

    printf " all good!\n"; 


*/
