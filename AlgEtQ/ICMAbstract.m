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

declare verbose ICMAbstract, 3;

declare attributes AlgEtQOrd: ICMAbstractRep;

declare type AlgEtQICM[AlgEtQICMElt];

declare attributes AlgEtQICM : Order,
                               ExtensionMapsPics, // an associative array with extension maps Pic(S)->Pic(T)
                               Map,
                               One;

declare attributes AlgEtQICMElt : Parent,
                                  WEClass,
                                  PicClass,
                                  IsOne,
                                  IsInvertibleInMultiplicatorRing;

///## Abstract representation of the ideal class monoid
/// The second way method to compute the ideal class monoid of an order $R$, returns an abstract representation of $\mathcal{I}(R)$ with type `AlgEtQICM` together with a map to a set of representatives. Each class has type `AlgEtQICMElt`, and it is internally represented as pair consisting of a weak equivalence class (of type `AlgEtQWECMElt`) and an element of the representation of $\mathrm{Pic}(S)$ (as an abstract abelian group), where $S$ is the corresponding multiplicator ring. This representation is more efficient than the previous one, since it does not need to compute in advance and store a representative for each ideal class.
/// Ideal classes can be created using the coercion operator `!` starting from an overorder $S$ of $R$ or a fractional $S$-ideal. Ideal classes can be multiplied using the operator `*`.

//////////////////
// AlgEtQICMElt //
//////////////////

///hide-all
function CreateAlgEtQICMElt(icm,yw,y_inv,pT)
// given the ideal class monoid of R, a weak equialence class of R, an element of the (abstract) picard gruop of an overorder of R, with its map pT:=Pic(T)-> ideals, it creates the corresponding ideal class.
    x:=New(AlgEtQICMElt);
    x`Parent:=icm;
    x`WEClass:=yw;
    x`PicClass:=<y_inv,pT>;
    return x;
end function;

intrinsic Print(x::AlgEtQICMElt)
{Print the element.}
    printf "ideal class over the order %o", Order(x);
end intrinsic;

intrinsic IsCoercible(icm::AlgEtQICM, x::.) -> BoolElt, .
{Return whether the element is coercible into icm and the result of the coercion if so.}
    if Parent(x) cmpeq icm then
        return true,x;
    elif Type(x) eq AlgEtQICMElt and Order(icm) subset Order(Parent(x)) then
        return true,(Map(Parent(x))(x))@@Map(icm); 
    elif Type(x) eq AlgEtQIdl and Order(icm) subset Order(x) then
        return true,(Order(icm)!!x)@@Map(icm); 
    elif Type(x) eq AlgEtQOrd and Order(icm) subset x then
        return true,(Order(icm)!!OneIdeal(x))@@Map(icm); 
    else 
        return false,"";
    end if;
end intrinsic;
///hide-none

/// Coerce $x$ into $\mathrm{icm}$, when possible.
intrinsic '!'(icm::AlgEtQICM, x::.) -> AlgEtQICMElt
{Coerce x into icm, when possible.}
    bool,x:=IsCoercible(icm,x);
    require bool : "The element cannot be coerced in the ideal class monoid.";
    return x;
end intrinsic;

/// Return whether $x$ is in $\mathrm{icm}$.
intrinsic 'in'(x::AlgEtQICMElt,icm::AlgEtQICM) -> BoolElt
{Returns whether x is in icm.}
    return Parent(x) cmpeq icm;
end intrinsic;

/// Returns the parent of $x$, that is, the ideal class monoid it belongs to.
intrinsic Parent(x::AlgEtQICMElt)->AlgEtQICM
{Returns the parent of x, that is, the ideal class monoid it belongs to.}
    return x`Parent;
end intrinsic;

/// Returns the weak equivalence class of x.
intrinsic WEClass(x::AlgEtQICMElt)->AlgEtQWECMElt
{Returns the weak equivalence class of x.}
    return x`WEClass;
end intrinsic;

/// Returns the elements of the abstract abelian group representing $\mathrm{Pic}(T)$, where $T$ is the multiplicator ring of x, together with the map from $\mathrm{Pic}(T)$ to the set of representing ideals.
intrinsic PicClass(x::AlgEtQICMElt)->GrpAbElt,Map
{Returns the elements of the PicardGroup Pic(T) of the multiplicator ring T of x together with the map from Pic(T) to the set of representing ideals.}
    return Explode(x`PicClass);
end intrinsic;

intrinsic Ideal(x::AlgEtQICMElt)->AlgEtQIdl
{Returns a deterministically computed ideal representing the given ideal class. This is not stored in an attribute to save memory.}
    m:=Map(Parent(x));
    return m(x);
end intrinsic;

intrinsic MultiplicatorRing(x::AlgEtQICMElt)->AlgEtQOrd
{Returns the multiplicator ring of the given ideal class.}
    return MultiplicatorRing(WEClass(x));
end intrinsic;

intrinsic 'eq'(x::AlgEtQICMElt,y::AlgEtQICMElt)->BoolElt
{Returns whether the two elements define the same class.}
    return WEClass(x) eq WEClass(y) and PicClass(x) eq PicClass(y);
end intrinsic;

//////////////////////////////////////////////
// Operations and extension map between Pic //
//////////////////////////////////////////////

ext_map_Pic:=function(icm,S)
// this function takes as input the icm of an order R, and an overorder S of R.
// it returns the extension map Pic(R)->Pic(S). 
// this populates the attribute ExtensionMapsPics of icm
    if not assigned icm`ExtensionMapsPics then
        icm`ExtensionMapsPics:=AssociativeArray();
    end if;
    test,map:=IsDefined(icm`ExtensionMapsPics,S);
    if not test then
        icm`ExtensionMapsPics[S]:=ExtensionHomPicardGroups(Order(icm),S);
    end if;
    return icm`ExtensionMapsPics[S];
end function;

intrinsic '*'(x::AlgEtQICMElt,y::AlgEtQICMElt)->AlgEtQICMElt
{Returns ideal class corresponding to the product of the two given ideal classes.}
    icm:=Parent(x);
    require icm eq Parent(y) : "The classes do not blong to the same ideal class monoid";
    Sx:=MultiplicatorRing(x);
    Sy:=MultiplicatorRing(y);
    wxy:=WEClass(x)*WEClass(y);
    Sxy:=MultiplicatorRing(wxy);
    _,pxy:=PicardGroup(Sxy);
    ext_x:=ext_map_Pic(icm,Sx);
    ext_y:=ext_map_Pic(icm,Sy);
    ext_xy:=ext_map_Pic(icm,Sxy);
    pic_xy:=PicClass(x)@@ext_x@ext_xy+PicClass(y)@@ext_y@ext_xy;
    return CreateAlgEtQICMElt(icm,wxy,pic_xy,pxy);
end intrinsic;
 
/// Given an ideal class $x$ and a non-negative integer $n$, returns the ideal class $x^n$.
intrinsic '^'(x::AlgEtQICMElt,n::RngIntElt)->AlgEtQICMElt
{Given an ideal class x and a non-negative integer n, returns the ideal class x^n.}
    require n ge 0 : "The integer must be non-negative.";
    if n eq 0 then
        return One(Parent(x));
    elif n eq 1 then
        return x;
    else
        icm:=Parent(x);
        S:=MultiplicatorRing(x);
        w_xn:=WEClass(x)^n;
        Sn:=MultiplicatorRing(w_xn);
        _,pn:=PicardGroup(Sn);
        ext_x:=ext_map_Pic(icm,S);
        ext_xn:=ext_map_Pic(icm,Sn);
        pic_xn:=n*(PicClass(x)@@ext_x@ext_xn); // additive notation for pic
        return CreateAlgEtQICMElt(icm,w_xn,pic_xn,pn);
    end if;
end intrinsic;

/////////
// icm //
/////////

/// Returns the ideal class monoid of the given order together with a map (with preimages) sending each class to a representative (determined by `WeakEquivalenceClassMonoid` and `PicardGroup`).
intrinsic IdealClassMonoidAbstract(R::AlgEtQOrd) -> AlgEtQICM,Map
{Returns the ideal class monoid of the given order together with a map (with preimages) sending each class to a representative (determined by WeakEquivalenceClassMonoid and PicardGroup).}
    if not assigned R`ICMAbstractRep then
        icm:=New(AlgEtQICM);
        W,w:=WeakEquivalenceClassMonoidAbstract(R);
        icm`Order:=R;

        img:=function(x)
            xw:=WEClass(x);
            xinv,pT:=PicClass(x);
            return w(xw)*(R!!pT(xinv));
        end function;
        preimg:=function(y)
            assert Order(y) eq R;
            yw:=y@@w;
            y_inv:=ColonIdeal(y,Ideal(yw)); // Ideal(yw)*y_inv = y
            T:=MultiplicatorRing(y);
            _,pT:=PicardGroup(T);
            y_inv:=(T!!y_inv)@@pT;
            return CreateAlgEtQICMElt(icm,yw,y_inv,pT);
        end function;
        m:=map<icm->PowerStructure(AlgEtQIdl)|x:->img(x) ,
                                              y:->preimg(y)>;
        icm`Map:=m;
        R`ICMAbstractRep:=icm;
    end if;
    return R`ICMAbstractRep,R`ICMAbstractRep`Map;
end intrinsic;

////////////////////////////////////
// AlgEtQICMElt : special classes //
////////////////////////////////////

intrinsic IsOne(x::AlgEtQICMElt)->BoolElt
{Returns whether the given ideal class is the neutral element of the ideal class monoid it belongs to.}
    if not assigned x`IsOne then
        x`IsOne:=x eq One(Parent(x));
    end if;
    return x`IsOne;
end intrinsic;

intrinsic IsInvertibleInMultiplicatorRing(x::AlgEtQICMElt)->BoolElt
{Returns whether given ideal class is invertible in its own multiplicator ring.}
    if not assigned x`IsInvertibleInMultiplicatorRing then
        x`IsInvertibleInMultiplicatorRing:=IsIdempotent(WEClass(x));
    end if;
    return x`IsInvertibleInMultiplicatorRing;
end intrinsic;

//////////////////////////////////////////////////
// AlgEtQICM : attributes and basic properties //
/////////////////////////////////////////////////

///hide-all
intrinsic Print(W::AlgEtQICM)
{Print the ideal class monoid.}
    printf "Ideal class monoid of %o", Order(W);
end intrinsic;
///hide-none

intrinsic Order(icm::AlgEtQICM)->AlgEtQOrd
{Returns the order of given ideal class monoid.}
    return icm`Order;
end intrinsic;

intrinsic Map(icm::AlgEtQICM)->Map
{Returns the map from the given ideal class monoid to the set of ideals which returns the representative of each class.}
    return icm`Map;
end intrinsic;

intrinsic 'eq'(icm1::AlgEtQICM,icm2::AlgEtQICM)->BoolElt
{Returns whether the two ideal class monoid are the same, that is, if the underlying orders are so.}
    return Order(icm1) eq Order(icm2);
end intrinsic;

intrinsic '#'(icm::AlgEtQICM)->RngInt
{Returns the size of the given ideal class monoid.}
    size:=0;
    W:=WeakEquivalenceClassMonoidAbstract(Order(icm));
    for T->v in Array(W) do
        size +:= #v*#PicardGroup(T);
    end for;
    return size;
end intrinsic;

intrinsic Classes(icm::AlgEtQICM)->SeqEnum[AlgEtQICMElt]
{Returns the sequence of the classes in the given ideal class monoid.}
    W:=WeakEquivalenceClassMonoidAbstract(Order(icm));
    output:=[];
    for T->WT in Array(W) do
        PT,pT:=PicardGroup(T);
        output cat:=[CreateAlgEtQICMElt(icm,w,z,pT):w in WT, z in PT];
    end for;
    return output;
end intrinsic;

intrinsic Representatives(icm::AlgEtQICM)->SeqEnum[AlgEtQIdl]
{Returns the sequence of representatives of the classes in the given ideal class monoid.}
    m:=Map(icm);
    return [m(c):c in Classes(icm)];
end intrinsic;

intrinsic One(icm::AlgEtQICM)->AlgEtQICMElt
{Returns the neutral element of the given ideal class monoid.}
    if not assigned icm`One then
        m:=Map(icm);
        icm`One:=OneIdeal(Order(icm))@@m;
    end if;
    return icm`One;
end intrinsic;

intrinsic Random(icm::AlgEtQICM)->AlgEtQICMElt
{Returns a random element of the given ideal class monoid.}
    return Random(Classes(icm));
end intrinsic;

///# Example 7 TODO

/* TESTS

    printf "### Testing ICMAbstract:";
	AttachSpec("~/AlgEt/spec");
    SetClassGroupBounds("GRH");

	_<x>:=PolynomialRing(Integers());
    f:=x^4+291*x^3-988*x^2-1000*x-1000;
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    icm,w:=IdealClassMonoidAbstract(E);
    assert forall{x:x in Classes(icm)| (w(x)*RandomUnit(K))@@w eq x};
    assert One(icm) eq icm!OneIdeal(E);
    assert IsOne(icm!(RandomUnit(K)*OneIdeal(E)));
    assert forall{T:T in OverOrders(E) | IsInvertibleInMultiplicatorRing(icm!T)};
    assert #IdealClassMonoidAbstract(MaximalOrder(K)) eq #PicardGroup(MaximalOrder(K));
    assert Random(icm) in icm;
    _:=Random(icm)*Random(icm);
    _:=[Random(icm)^i:i in [1..100]];

	_<x>:=PolynomialRing(Integers());
    f:=x^3+31*x^2+43*x+77;
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    icm,w:=IdealClassMonoidAbstract(E);
    assert forall{x:x in Classes(icm)| (w(x)*RandomUnit(K))@@w eq x};
    assert One(icm) eq icm!OneIdeal(E);
    assert IsOne(icm!(RandomUnit(K)*OneIdeal(E)));
    assert forall{T:T in OverOrders(E) | IsInvertibleInMultiplicatorRing(icm!T)};
    assert #IdealClassMonoidAbstract(MaximalOrder(K)) eq #PicardGroup(MaximalOrder(K));
    assert Random(icm) in icm;
    _:=Random(icm)*Random(icm);
    _:=[Random(icm)^i:i in [1..100]];
    _:=[x*y:x,y in Classes(icm)];

    L:=RationalsAsNumberField();
    K:=EtaleAlgebra([L,L,L]);
    a:=PrimitiveElement(K);
    E:=Order([a]);
    icm:=IdealClassMonoidAbstract(E);
    O:=MaximalOrder(K);
    assert #IdealClassMonoidAbstract(O) eq 1;
    assert not IsMaximal(E);
    assert One(icm) eq icm!OneIdeal(E);
    assert IsOne(icm!(RandomUnit(K)*OneIdeal(E)));
    _:=Random(icm)*Random(icm);
    _:=[Random(icm)^i:i in [1..100]];
    _:=[x*y:x,y in Classes(icm)];

    SetAssertions(1);
    printf " all good!"; 
*/

