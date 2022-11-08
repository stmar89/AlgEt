/* vim: set syntax=magma :*/

freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

//declare verbose ?????, 1;

/*TODO:

*/

RemoveBlanks:=function(str)
// given a string str removes the blank spaces
    return &cat(Split(str," "));
end function;

//------------
// Print to string
//------------

intrinsic PrintSeqAlgEtElt(seq::SeqEnum[AlgEtElt]) -> SeqEnum,MonStgElt
{ Given a sequence of elements of an AlgEt, returns a sequence of tuples of sequence of integers that can be coerced into the original algebra to obtain the input sequece. As a second output it reuturns a string that can be printed to file.}
    seq:=[ < Eltseq(c) : c in Components(elt) > : elt in seq ];
    str:=RemoveBlanks(Sprint(seq));
    return seq,str;
end intrinsic;

intrinsic PrintWKICM(R::AlgEtOrd) -> MonStgElt
{ Given an order R in an AlgEt, it returns a string that contains the weak equivalence classes of R, sorted by multiplicator ring. In particular, the overorders of R can be recovered fro this string. Such string can be easily printed to file. To load the string, after using Read() on the file, use the intrinsic LoadWKICM }
    str:="<\n";
    A:=Algebra(R);
    nf:=NumberFields(A);
    nf_poly:=[ Coefficients((DefiningPolynomial(K))) : K in nf ];
    str cat:=RemoveBlanks(Sprint(nf_poly)) cat ",\n";
    oo:=FindOverOrders(R);
    // we make sure that R is oo[1] and the max order O is in the last position
    iR:=Index(oo,R);
    if iR ne 1 then
        temp:=oo[1];
        oo[1]:=R;
        oo[iR]:=temp;
    end if;
    O:=MaximalOrder(A);
    iO:=Index(oo,O);
    if iO ne #oo then
        temp:=oo[#oo];
        oo[#oo]:=O;
        oo[iO]:=temp;
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
            _,strI:=PrintSeqAlgEtElt(ZBasis(I));
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

intrinsic LoadWKICM(str::MonStgElt) -> AlgEtOrd
{ Given a string produced with PrintWKICM it returns the corresponding order R. In the attributes of R, its algebra and its overorders one can find the weak equivalence classes. This can be recovered with the approriate intrinsics. }

    data:=eval(str);
    PP:=PolynomialRing(Rationals());
    ff:=[ PP!f : f in data[1]];
    A:=EtaleAlgebra([NumberField(f) : f in ff ]);
    wk:=data[2];
    O:=Order([ A ! s : s in wk[#wk][1]]); // the first one is maximal
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
        R:=Order([ A! s : s in wk[1][1]]); // the first one is R
    else 
        R:=O; //to save attributes
    end if;
    ooR:={@ @};
    wkR:=[];
    for iS->dataS in wk do
        if iS eq 1 then
            S:=R;
        elif iS eq #wk then
            S:=O;
        else
            S:=Order([ A ! s : s in dataS[1]]);
        end if;
        S`IsGorenstein:=#dataS eq 1;
        Include(~ooR,S);
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



/* TEST
    
    AttachSpec("~/packages_github/AlgEt/spec");
    _<x>:=PolynomialRing(Integers());
    f:=x^6 - 3*x^5 - 3*x^4 + 65*x^3 - 48*x^2 - 768*x + 4096;
    A:=EtaleAlgebra(f);
    E:=EquationOrder(A);
    seq,str:=PrintSeqAlgEtElt(ZBasis(E));
    assert Order([ A! s : s in eval(str)]) eq E;

    
    AttachSpec("~/packages_github/AlgEt/spec");
    _<x>:=PolynomialRing(Integers());
    f:=x^8+16;
    A:=EtaleAlgebra(f);
    O:=MaximalOrder(A);
    time str:=PrintWKICM(O);
    time O1:=LoadWKICM(str);

    AttachSpec("~/packages_github/AlgEt/spec");
    _<x>:=PolynomialRing(Integers());
    f:=x^8+16;
    A:=EtaleAlgebra(f);
    F:=PrimitiveElement(A);
    R:=Order([F,2/F]);
    time str:=PrintWKICM(R);
    time R1:=LoadWKICM(str);
    assert #WKICM(R) eq #WKICM(R1);
    assert #FindOverOrders(R) eq #FindOverOrders(R1);

*/
