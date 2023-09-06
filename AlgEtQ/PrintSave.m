/* vim: set syntax=magma :*/

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
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

RemoveBlanks:=function(str)
// given a string str removes the blank spaces
    return &cat(Split(str," "));
end function;

//------------
// Print to string
//------------

intrinsic PrintSeqAlgEtQElt(seq::SeqEnum[AlgEtQElt]) -> SeqEnum,MonStgElt
{Given a sequence of elements of an AlgEtQ, returns a sequence of tuples of sequence of integers that can be coerced into the original algebra to obtain the input sequece. As a second output it returns a string that can be printed to file.}
    seq:=[ < Eltseq(c) : c in Components(elt) > : elt in seq ];
    str:=RemoveBlanks(Sprint(seq));
    return seq,str;
end intrinsic;

intrinsic PrintWKICM(R::AlgEtQOrd) -> MonStgElt
{Given an order R in an AlgEtQ, it returns a string that contains the weak equivalence classes of R, sorted by multiplicator ring. In particular, the overorders of R can be recovered from this string. Such string can be easily printed to file. To load the string, after using Read() on the file, use the intrinsic LoadWKICM.}
    str:="<\n";
    A:=Algebra(R);
    nf:=Components(A);
    nf_poly:=[ Coefficients((DefiningPolynomial(K))) : K in nf ];
    str cat:=RemoveBlanks(Sprint(nf_poly)) cat ",\n";
    oo:=FindOverOrders(R);
    // we make sure that R is oo[1]
    if R ne oo[1] then
        ind:=Index(oo,R);
        assert ind ne 0;
        S:=oo[1];
        oo[ind]:=S;
        oo[1]:=R;
    end if;
    // init the string
    str cat:="<\n";
    for iS->S in oo do
        wkS:=WKICM_bar(S);
        // we sort it in such a way that the invertible class is the first and represented by 1*S
        // this is not done in WkClasses.m
        invcl:=[ IsInvertible(I) select 1 else 0 : I in wkS ];
        assert &+invcl eq 1; //only one inv class
        ind:=Index(invcl,1);
        wkS:=[ OneIdeal(S) ] cat Remove(wkS,ind);
        // produce the string
        str cat:="<\n";
        for iI->I in wkS do
            _,strI:=PrintSeqAlgEtQElt(ZBasis(I));
            if iI ne #wkS then
                str cat:=strI cat ",\n";
            else
                str cat:=strI cat "\n";
            end if;
        end for;
        if iS ne #oo then
            str cat:= ">,\n";
        else
            str cat:= ">\n";
        end if;
    end for;
    str cat:= ">\n>";
    return str;
end intrinsic;

intrinsic LoadWKICM(str::MonStgElt) -> AlgEtQOrd
{Given a string produced with PrintWKICM it returns the corresponding order R. In the attributes of R, its algebra and its overorders one can find the weak equivalence classes. These can be recovered with the approriate intrinsics.}

    data:=eval(str);
    PP:=PolynomialRing(Rationals());
    ff:=[ PP!f : f in data[1]];
    A:=EtaleAlgebra([NumberField(f) : f in ff ]);
    wk:=data[2];

    ooR:=[ Order([ A ! s : s in wk[j][1]]) : j in [1..#wk] ];
    indices:=[ Index(S,ooR[1]) : S in ooR ];
    max,indO:=Max(indices);
    assert #[ i : i in indices | i eq max ] eq 1; //sanity check
    O:=ooR[indO]; // this is the maximal order of A
    O`IsMaximal:=true;
    test,OasProd:=IsProductOfOrders(O);
    assert test;
    for i in [1..#OasProd] do
        OL:=OasProd[i];
        OL`Maximal:=true;
        OL`MaximalOrder:=OL;
    end for;
    O`IsGorenstein:=true;
    A`MaximalOrder:=O;
    if #wk gt 1 then
        R:=ooR[1]; // the first one is R
    else 
        R:=O; //to save attributes
        assert indO eq 1; //sanity check
    end if;
    wkR:=[];
    for iS->dataS in wk do
        S:=ooR[iS];
        S`IsGorenstein:=#dataS eq 1; //only one weak equivalence class
        wkS:=[];
        for iI->I in dataS do
            zbI:=[A!s : s in I];
            IS:=Ideal(S,zbI);
            IS`MultiplicatorRing:=S;
            IS`ZBasis:=zbI;
            IS`IsInvertible:=iI eq 1; //only the first is invertible
            Append(~wkS,IS);
            IR:=R!!IS;
            Append(~wkR,IR); //attributes are already moved
        end for;
        S`WKICM_bar:=wkS;
    end for;
    R`OverOrders:=ooR;
    R`WKICM:=wkR;
    return R;
end intrinsic;



/* TESTS
    
    printf "### Testing Print Saving:";
    //AttachSpec("~/packages_github/AlgEt/spec");
    _<x>:=PolynomialRing(Integers());
    f:=x^6 - 3*x^5 - 3*x^4 + 65*x^3 - 48*x^2 - 768*x + 4096;
    A:=EtaleAlgebra(f);
    E:=EquationOrder(A);
    seq,str:=PrintSeqAlgEtQElt(ZBasis(E));
    assert Order([ A! s : s in eval(str)]) eq E;
    printf ".";

    //AttachSpec("~/packages_github/AlgEt/spec");
    _<x>:=PolynomialRing(Integers());
    f:=x^8+16;
    A:=EtaleAlgebra(f);
    O:=MaximalOrder(A);
    str:=PrintWKICM(O);
    O1:=LoadWKICM(str);
    printf ".";

    //AttachSpec("~/packages_github/AlgEt/spec");
    _<x>:=PolynomialRing(Integers());
    f:=x^8+16;
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    R:=Order([F,2/F]);
    str:=PrintWKICM(R);
    R1:=LoadWKICM(str);
    assert #WKICM(R) eq #WKICM(R1);
    assert #FindOverOrders(R) eq #FindOverOrders(R1);
    printf ".";

    printf " all good!\n"; 

*/
