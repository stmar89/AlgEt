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
// Copyright 2023, S. Marseglia
/////////////////////////////////////////////////////

freeze;

declare verbose WkAbstract, 3;

declare attributes AlgEtQOrd: WKICMAbstractRep;

declare type AlgEtQWECM[AlgEtQWECMElt];

declare attributes AlgEtQWECM : Order,
                                Array,
                                Map,
                                One,
                                MultiplicationTable;

declare attributes AlgEtQWECMElt : Parent,
                                   Ideal,
                                   IsOne,
                                   IsIdempotent;

///////////////////////////////////
// Creation of elements and wkcm //
///////////////////////////////////

function CreateAlgEtQWECMElt(W,I)
// a the weak equiv class monoid of an order R, of type AlgEtQWECM and a fractional R-ideal, it creates an element of the monoid
// to be used only when creating the weak equivalence class monoid
    x:=New(AlgEtQWECMElt);
    x`Parent:=W;
    x`Ideal:=I;
    return x;
end function;

intrinsic WeakEquivalenceClassMonoidAbstract(R::AlgEtQOrd : Method:="Auto") -> AlgEtQWECM,Map
{Returns the weak equivalence class monoid W of R together with a map w (with preimages) sending each class of W to a representative (determined by WeakEquivalenceClassMonoid).}
    if not assigned R`WKICMAbstractRep then
        if (not assigned R`OverOrders) or exists{T:T in R`OverOrders|not assigned T`WKICM_bar} then
            // this populates everything faster than computing separately the overorders and the wk_icm_bar's 
            _:=WeakEquivalenceClassMonoid(R);
        end if;
        oo:=OverOrders(R);
        W:=New(AlgEtQWECM);
        W`Order:=R;

        arr:=AssociativeArray();
        for T in oo do 
            arr[T]:=[CreateAlgEtQWECMElt(W,R!!I) : I in WKICM_bar(T)];
        end for;
        W`Array:=arr;
        // we now defined the map
        dlp:=function(y)
            assert Order(y) eq R ;
            T:=MultiplicatorRing(y);
            assert exists(I){w:w in arr[T]|IsWeaklyEquivalent(Ideal(w),y)};
            return I;
        end function;

        w:=map<W->PowerStructure(AlgEtQIdl)|x:->Ideal(x),
                                            y:->dlp(y)>;
        W`Map:=w;
        R`WKICMAbstractRep:=W;
    end if;
    return R`WKICMAbstractRep,R`WKICMAbstractRep`Map;
end intrinsic;

///////////////////
// AlgEtQWECMElt //
///////////////////

intrinsic Print(x::AlgEtQWECMElt)
{Print the element.}
    printf "weak equivalence class of the ideal %o", Ideal(x);
end intrinsic;

intrinsic IsCoercible(W::AlgEtQWECM, x::.) -> BoolElt, .
{Return whether the element is coercible into W and the result of the coercion if so.}
    if Parent(x) cmpeq W then
        return true,x;
    elif Type(x) eq AlgEtQWECMElt and Order(W) subset Order(Parent(x)) then
        return true,(Order(W)!!Ideal(x))@@Map(W); 
    elif Type(x) eq AlgEtQIdl and Order(W) subset Order(x) then
        return true,(Order(W)!!x)@@Map(W); 
    elif Type(x) eq AlgEtQOrd and Order(W) subset x then
        return true,(Order(W)!!OneIdeal(x))@@Map(W); 
    else 
        return false,"";
    end if;
end intrinsic;

intrinsic '!'(W::AlgEtQWECM, x::.) -> AlgEtQWECMElt
{Coerce x into W.}
    bool,x:=IsCoercible(W,x);
    require bool : "The element cannot be coerced in the weak equivalence class monoid.";
    return x;
end intrinsic;

intrinsic 'in'(x::AlgEtQWECMElt,W::AlgEtQWECM) -> BoolElt
{Returns whether x is in W.}
    return Parent(x) cmpeq W;
end intrinsic;

intrinsic Parent(x::AlgEtQWECMElt)->AlgEtQWECM
{Returns the parent of x, that is, the weak equivalence class monoid it belongs to.}
    return x`Parent;
end intrinsic;

intrinsic Ideal(x::AlgEtQWECMElt)->AlgEtQIdl
{Returns the ideal representing x.}
    return x`Ideal;
end intrinsic;

intrinsic MultiplicatorRing(x::AlgEtQWECMElt)->AlgEtQOrd
{Returns the multiplicator ring of x.}
    return MultiplicatorRing(Ideal(x));
end intrinsic;

intrinsic 'eq'(x::AlgEtQWECMElt,y::AlgEtQWECMElt)->BoolElt
{Returns whether the two elements define the same class.}
    I:=Ideal(x);
    J:=Ideal(y);
    if I eq J then
        return true;
    else
        assert3 not IsWeaklyEquivalent(I,J);
        return false;
    end if;
end intrinsic;

intrinsic SetRepresentative(x::AlgEtQWECMElt,J::AlgEtQIdl)
{It changes the Ideal attribute of x which consisits of the ideal representing the weak equivalence class to J, which must be weakly equivalent to the one of x.}
    require IsWeaklyEquivalent(Ideal(x),J) : "The given ideal does not define the same class";
    x`Ideal:=J;
end intrinsic;

////////////////////////////////////////
// Operation and Multiplication Table //
///////////////////////////////////////

intrinsic MultiplicationTable(W::AlgEtQWECM)->Assoc
{Computes the multiplication table of W.}
    if not assigned W`MultiplicationTable then
        W`MultiplicationTable:=AssociativeArray();
    end if;
    kk:=Keys(W`MultiplicationTable);
    cl:=Classes(W);
    for x,y in cl do
        //this fills all the missing entries. No operation is done for the already computed ones.
        if {x,y} notin kk then
            _:=x*y;
        end if;
    end for;
    return W`MultiplicationTable;
end intrinsic;

intrinsic '*'(x::AlgEtQWECMElt,y::AlgEtQWECMElt)->AlgEtQWECMElt
{Returns weak equivalence class corresponding to the product.}
    W:=Parent(x);
    require W eq Parent(y) : "The classes do not blong to the same weak equivalence class monoid";
    if not assigned W`MultiplicationTable then
        W`MultiplicationTable:=AssociativeArray();
    end if;
    test,z:= IsDefined(W`MultiplicationTable,{x,y});
    if test then
        return z;
    else
        if IsOne(x) then
            z:=y;
        elif IsOne(y) then
            z:=x;
        else
            z:=W!(Ideal(x)*Ideal(y));
        end if;
        W`MultiplicationTable[{x,y}]:=z;
    end if;
    return z;
end intrinsic;
 
intrinsic '^'(x::AlgEtQWECMElt,n::RngIntElt)->AlgEtQWECMElt
{Given a weak equivalence class x and a non-negative integer n, returns the weak equivalence class x^n.}
    require n ge 0 : "The integer must be non-negative.";
    if n eq 0 then
        return One(Parent(x));
    elif n eq 1 then
        return x;
    elif assigned x`IsIdempotent and IsIdempotent(x) then
        return x;
    elif n eq 2 then
        return x*x; //using the multiplication table
    else
        return Parent(x)!(Ideal(x)^n);
    end if;
end intrinsic;


/////////////////////////////////////
// AlgEtQWECMElt : special classes //
////////////////////////////////////

intrinsic IsOne(x::AlgEtQWECMElt)->BoolElt
{Returns whether x is the neutral element of the weak equivalence class monoid it belongs to.}
    return x eq One(Parent(x));
end intrinsic;

intrinsic IsIdempotent(x::AlgEtQWECMElt)->BoolElt
{Returns whether x is an idempotent of weak equivalence class monoid it belongs to.}
    if not assigned x`IsIdempotent then
        x`IsIdempotent:=x*x eq x;
    end if;
    return x`IsIdempotent;
end intrinsic;

//////////////////////////////////////////////////
// AlgEtQWECM : attributes and basic properties //
/////////////////////////////////////////////////

intrinsic Print(W::AlgEtQWECM)
{Print the weak equivalence class monoid.}
    printf "Weak equivalence class monoid of %o", Order(W);
end intrinsic;

intrinsic Order(W::AlgEtQWECM)->AlgEtQOrd
{Returns the order of W.}
    return W`Order;
end intrinsic;

intrinsic Array(W::AlgEtQWECM)->Assoc
{Return the underlying associative array of W, which is indexed by the overorder of R and values given by the weak equivalence classes with prescribed multiplicator ring.}
    return W`Array;
end intrinsic;

intrinsic Map(W::AlgEtQWECM)->Map
{Returns the map from W to the set of ideals which returns the representative of each class.}
    return W`Map;
end intrinsic;

intrinsic 'eq'(W::AlgEtQWECM,WW::AlgEtQWECM)->BoolElt
{Returns whether the two weak equivalence class monoid are the same, that is, if the underlying orders are.}
    return Order(W) eq Order(WW);
end intrinsic;

intrinsic '#'(W::AlgEtQWECM)->RngInt
{Returns the size of W.}
    size:=0;
    for k->v in Array(W) do
        size +:= #v;
    end for;
    return size;
end intrinsic;

intrinsic Classes(W::AlgEtQWECM)->SeqEnum[AlgEtQWECMElt]
{Returns the sequence of the classes in W.}
    arr:=Array(W);
    return &cat[ arr[T] : T in Keys(arr)];
end intrinsic;

intrinsic Representatives(W::AlgEtQWECM)->SeqEnum[AlgEtQIdl]
{Returns the sequence of representatives of the classes in W.}
    w:=Map(W);
    return [w(c):c in Classes(W)];
end intrinsic;

intrinsic One(W::AlgEtQWECM)->AlgEtQWECMElt
{Returns the neutral element of W.}
    if not assigned W`One then
        w:=Map(W);
        W`One:=OneIdeal(Order(W))@@w;
    end if;
    return W`One;
end intrinsic;

intrinsic Random(W::AlgEtQWECM)->AlgEtQWECMElt
{Returns a random element of W.}
    return Random(Classes(W));
end intrinsic;

intrinsic Idempotents(W::AlgEtQWECM)->SeqEnum[AlgEtQWECMElt]
{Returns the sequence of the classes in W which are idempotent.}
    return [x:x in Classes(W)|IsIdempotent(x)];
end intrinsic;

////////////////////////////////
// AlgEtQWECM : localizations //
////////////////////////////////

intrinsic Localization(W::AlgEtQWECM,P::AlgEtQIdl)->AlgEtQWECM
{Given a weak equivalence class monoid of an order R in the etale algebra K and a maximal ideal P of R, returns the weak equivalance class monoid of the unique overorder of R which is locally equal to R at P and locally maxiaml at every other maximal ideal. This order is R+P^kO, where O is the maximal order of K and k is a non-negative integer big enough.}
    R:=Order(W);
    require Order(P) eq R and IsPrime(P) : "P needs to be a maximal ideal of the underlying order of W";
    O:=MaximalOrder(Algebra(R));
    _,p:=IsPrimePower(Index(R,P));
    k:=Valuation(Index(O,R),p);
    RP:=Order(ZBasis(R) cat ZBasis(O!!P^k));
    return WeakEquivalenceClassMonoidAbstract(RP);
end intrinsic;

//////////////////////////////////
// AlgEtQWECM : generating sets //
/////////////////////////////////

intrinsic IsGeneratingSet(W::AlgEtQWECM,seq::SeqEnum[AlgEtQWECMElt])->BoolElt
{Given the weak equivalence class monoid of R and a sequence of classes in W, returns whether the sequence is a generating set of W.}
    cl:=Seqset(seq);
    require forall{w:w in seq|w in W} : "The sequence does not consists of elements of W.";
    submonoid_old:=cl;
    // TODO I am redoing the same multiplication over and over
    repeat
        submonoid_new:=submonoid_old;
        N:=#submonoid_old;
        for x,y in submonoid_old do
           Include(~submonoid_new,x*y); 
        end for;
        submonoid_old:=submonoid_new;
    until #submonoid_old eq N;
    return N eq #W;
end intrinsic;

intrinsic IsGeneratingSet(W::AlgEtQWECM,seq::SeqEnum[AlgEtQIdl])->BoolElt
{Given the weak equivalence class monoid of R and a sequence of fractional R-ideals, returns whether the sequence represents a generating set of W.}
    w:=Map(W);
    return IsGeneratingSet(W,[I@@w:I in Seqset(seq)]);
end intrinsic;

/* TESTS
TODO Test on QxQ or KxK

    printf "### Testing WkAbstract:";
	AttachSpec("~/AlgEt/spec");
	SetAssertions(2);
    SetClassGroupBounds("GRH");
	_<x>:=PolynomialRing(Integers());
    f:=x^4-1000*x^3-1000*x^2-1000*x-1000;
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    W,w:=WeakEquivalenceClassMonoidAbstract(E);
    x:=Random(W);
    xx:=x;
    SetRepresentative(x,Random(K)*Ideal(x));
    assert x eq xx;
    wk:=WKICM(E);
    wk1:=[w*Random(K):w in wk];
    assert #W eq 25;
    assert Seqset(Classes(W)) eq {x@@w : x in wk1};
    assert One(W) eq W!OneIdeal(E);
    assert IsOne(W!(Random(K)*OneIdeal(E)));
    assert IsIdempotent(W!MaximalOrder(K));
    assert #WeakEquivalenceClassMonoidAbstract(MaximalOrder(K)) eq 1;
    assert Random(W) in W;
    assert Seqset(Idempotents(W)) eq {W!T:T in OverOrders(E)};
    assert #W eq &*[#Localization(W,P) : P in SingularPrimes(E)];
    assert WeakEquivalenceClassMonoidAbstract(MaximalOrder(K)) eq Localization(W,PrimesAbove(11*E)[1]); // 11 is regular
    _:=Random(W)*Random(W);
    assert assigned W`MultiplicationTable;
    assert not IsGeneratingSet(W,Idempotents(W));
    assert Max([CohenMacaulayType(T):T in OverOrders(E)]) eq 2; //this implies the following line
    assert IsGeneratingSet(W,Idempotents(W) cat [W!TraceDualIdeal(T):T in OverOrders(E)]);
    _:=MultiplicationTable(W);

	_<x>:=PolynomialRing(Integers());
    f:=x^4+291*x^3-988*x^2-1000*x-1000;
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    W,w:=WeakEquivalenceClassMonoidAbstract(E);
    wk:=WKICM(E);
    P,p:=PicardGroup(E);
    inv:=[p(x):x in P];
    wk1:=[w*Random(inv):w in wk];
    assert #W eq 20;
    assert Seqset(Classes(W)) eq {x@@w : x in wk1};
    assert One(W) eq W!OneIdeal(E);
    assert IsOne(W!(Random(K)*OneIdeal(E)));
    assert IsIdempotent(W!MaximalOrder(K));
    assert #WeakEquivalenceClassMonoidAbstract(MaximalOrder(K)) eq 1;
    assert Random(W) in W;
    assert Seqset(Idempotents(W)) eq {W!T:T in OverOrders(E)};
    assert #W eq &*[#Localization(W,P) : P in SingularPrimes(E)];
    assert WeakEquivalenceClassMonoidAbstract(MaximalOrder(K)) eq Localization(W,PrimesAbove(11*E)[1]); // 11 is regular
    assert not IsGeneratingSet(W,Idempotents(W));
    assert Max([CohenMacaulayType(T):T in OverOrders(E)]) eq 2; //this implies the following line
    assert IsGeneratingSet(W,Idempotents(W) cat [W!TraceDualIdeal(T):T in OverOrders(E)]);

	_<x>:=PolynomialRing(Integers());
    f:=x^3+31*x^2+43*x+77;
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    W,w:=WeakEquivalenceClassMonoidAbstract(E);
    wk:=WKICM(E);
    P,p:=PicardGroup(E);
    inv:=[p(x):x in P];
    wk1:=[w*Random(inv):w in wk];
    assert #W eq 23;
    assert Seqset(Classes(W)) eq {x@@w : x in wk1};
    assert One(W) eq W!OneIdeal(E);
    assert IsOne(W!(Random(K)*OneIdeal(E)));
    assert IsIdempotent(W!MaximalOrder(K));
    assert #WeakEquivalenceClassMonoidAbstract(MaximalOrder(K)) eq 1;
    assert Random(W) in W;
    assert Seqset(Idempotents(W)) eq {W!T:T in OverOrders(E)};
    assert #W eq &*[#Localization(W,P) : P in SingularPrimes(E)];
    assert WeakEquivalenceClassMonoidAbstract(MaximalOrder(K)) eq Localization(W,PrimesAbove(11*E)[1]); // 11 is regular
    assert not IsGeneratingSet(W,Idempotents(W));
    assert Max([CohenMacaulayType(T):T in OverOrders(E)]) eq 2; //this implies the following line
    assert IsGeneratingSet(W,Idempotents(W) cat [W!TraceDualIdeal(T):T in OverOrders(E)]);

    SetAssertions(1);
    printf " all good!"; 
*/

