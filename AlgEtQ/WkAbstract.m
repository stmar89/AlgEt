/////////////////////////////////////////////////////
// Copyright 2025.
// Stefano Marseglia, stefano.marseglia89@gmail.com
// https://stmar89.github.io/index.html
// 
// Distributed under the terms of the CC-BY 4.0 licence.
// https://creativecommons.org/licenses/by/4.0/
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

///## Abstract representation of the weak equivalence class monoid.
/// The second method to compute $\mathcal{W}(R)$ returns an abstract representation of the monoid together with a map giving representatives. The monoid is of type `AlgEtQWECM` and each class is of type `AlgEtQWECMElt`. Classes can be created by feeding to the operator `!` a fractional $S$-ideal, when $S$ is an overorder of $R$ or an overorder of $R$. 
/// Whenever a class is created a representative is stored in an attribute.
/// This representative can be changed using the intrinsic `SetRepresentative`.
/// Classes can be multiplied using the operator `*`. Whenever two classes are multiplied, the result is stored in the attribute `MultiplicationTable` of the monoid.


///////////////////
// AlgEtQWECMElt //
///////////////////

function CreateAlgEtQWECMElt(W,I)
// a the weak equiv class monoid of an order R, of type AlgEtQWECM and a fractional R-ideal, it creates an element of the monoid
// to be used only when creating the weak equivalence class monoid
    x:=New(AlgEtQWECMElt);
    x`Parent:=W;
    x`Ideal:=I;
    return x;
end function;

///hide-all
intrinsic Print(x::AlgEtQWECMElt)
{Print the element.}
    printf "weak equivalence class of the ideal %o", Ideal(x);
end intrinsic;

intrinsic IsCoercible(W::AlgEtQWECM, x::.) -> BoolElt, .
{Returns whether the element is coercible into W and the result of the coercion if so.}
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
///hide-none

/// Given a weak equivalence class monoid $W$, coerces $x$ into $W$, when possible.
intrinsic '!'(W::AlgEtQWECM, x::.) -> AlgEtQWECMElt
{Given a weak equivalence class monoid W, coerces x into W, when possible.}
    bool,x:=IsCoercible(W,x);
    require bool : "The element cannot be coerced in the weak equivalence class monoid.";
    return x;
end intrinsic;

/// Returns whether the class $x$ is in the weak equivalence class monoid $W$.
intrinsic 'in'(x::AlgEtQWECMElt,W::AlgEtQWECM) -> BoolElt
{Returns whether the class x is in the weak equivalence class monoid W.}
    return Parent(x) cmpeq W;
end intrinsic;

/// Returns the weak equivalence class monoid to which the class $x$ belongs to.
intrinsic Parent(x::AlgEtQWECMElt)->AlgEtQWECM
{Returns the parent of x, that is, the weak equivalence class monoid it belongs to.}
    return x`Parent;
end intrinsic;

intrinsic Ideal(x::AlgEtQWECMElt)->AlgEtQIdl
{Returns the stored ideal representing the class.}
    return x`Ideal;
end intrinsic;

intrinsic MultiplicatorRing(x::AlgEtQWECMElt)->AlgEtQOrd
{Returns the multiplicator ring of the class.}
    return MultiplicatorRing(Ideal(x));
end intrinsic;

intrinsic 'eq'(x::AlgEtQWECMElt,y::AlgEtQWECMElt)->BoolElt
{Equality testing for classes.}
    I:=Ideal(x);
    J:=Ideal(y);
    if I eq J then
        return true;
    else
        assert3 not IsWeaklyEquivalent(I,J);
        return false;
    end if;
end intrinsic;

/// Given a class $x$ and a fractional ideal $J$ weakly equivalent to the stored representative, changes the attribute `Ideal` of the class to $J$.
intrinsic SetRepresentative(x::AlgEtQWECMElt,J::AlgEtQIdl)
{Given a class x and a fractional ideal J weakly equivalent to the stored representative, changes the attribute Ideal of the class to J.}
    require IsWeaklyEquivalent(Ideal(x),J) : "The given ideal does not define the same class";
    x`Ideal:=J;
end intrinsic;

////////////////////////////////////////
// Operation and Multiplication Table //
///////////////////////////////////////

/// Computes the multiplication table of $W$. It is returned as an associative array where, given two classes $x$ and $y$ of $W$, the value at the key $\{x,y\}$ is $x*y$.
intrinsic MultiplicationTable(W::AlgEtQWECM)->Assoc
{Computes the multiplication table of W. It is returned as an associative array where, given two classes x and y of W, the value at the key \{x,y\} is x*y.}
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
{Returns weak equivalence class corresponding to the product of the two given classes.}
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
 
/// Given a weak equivalence class $x$ and a non-negative integer $n$, returns the weak equivalence class $x^n$.
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
{Returns whether the given class is the neutral element of the weak equivalence class monoid it belongs to.}
    if not assigned x`IsOne then
        x`IsOne:=x eq One(Parent(x));
    end if;
    return x`IsOne;
end intrinsic;

intrinsic IsIdempotent(x::AlgEtQWECMElt)->BoolElt
{Returns whether the given class is an idempotent of weak equivalence class monoid it belongs to.}
    if not assigned x`IsIdempotent then
        x`IsIdempotent:=x*x eq x;
    end if;
    return x`IsIdempotent;
end intrinsic;

///////////////////////////////////////
// AlgEtQWECM : computation/creation //
///////////////////////////////////////

/// Given an order $R$, returns the weak equivalence class monoid $W$ of $R$ together with a map (with preimages) sending each class of $W$ to a representative (determined by WeakEquivalenceClassMonoid).}
intrinsic WeakEquivalenceClassMonoidAbstract(R::AlgEtQOrdd) -> AlgEtQWECM,Map
{Given an order R, returns the weak equivalence class monoid W of R together with a map (with preimages) sending each class of W to a representative (determined by WeakEquivalenceClassMonoid).}
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
        preimg:=function(y)
            assert Order(y) eq R ;
            T:=MultiplicatorRing(y);
            assert exists(I){w:w in arr[T]|IsWeaklyEquivalent(Ideal(w),y)};
            return I;
        end function;

        w:=map<W->PowerStructure(AlgEtQIdl)|x:->Ideal(x),
                                            y:->preimg(y)>;
        W`Map:=w;
        R`WKICMAbstractRep:=W;
    end if;
    return R`WKICMAbstractRep,R`WKICMAbstractRep`Map;
end intrinsic;

//////////////////////////////////////////////////
// AlgEtQWECM : attributes and basic properties //
/////////////////////////////////////////////////

///hide-all
intrinsic Print(W::AlgEtQWECM)
{Print the weak equivalence class monoid.}
    printf "Weak equivalence class monoid of %o", Order(W);
end intrinsic;
///hide-none

/// Returns the order of $W$.
intrinsic Order(W::AlgEtQWECM)->AlgEtQOrd
{Returns the order of W.}
    return W`Order;
end intrinsic;

/// Returns the underlying associative array of $W$, which is indexed by the overorders of $R$ and has values given by the weak equivalence classes with prescribed multiplicator ring.
intrinsic Array(W::AlgEtQWECM)->Assoc
{Returns the underlying associative array of W, which is indexed by the overorders of R and has values given by the weak equivalence classes with prescribed multiplicator ring.}
    return W`Array;
end intrinsic;

/// Returns the map from $W$ to the set of ideals which returns the representative of each class.
intrinsic Map(W::AlgEtQWECM)->Map
{Returns the map from W to the set of ideals which returns the representative of each class.}
    return W`Map;
end intrinsic;

intrinsic 'eq'(W::AlgEtQWECM,WW::AlgEtQWECM)->BoolElt
{Equality testing for weak equivalence class monoids, that is, if the underlying orders are equal.}
    return Order(W) eq Order(WW);
end intrinsic;

/// Returns the size of $W$.
intrinsic '#'(W::AlgEtQWECM)->RngInt
{Returns the size of W.}
    size:=0;
    for k->v in Array(W) do
        size +:= #v;
    end for;
    return size;
end intrinsic;

/// Returns the sequence of classes in $W$.
intrinsic Classes(W::AlgEtQWECM)->SeqEnum[AlgEtQWECMElt]
{Returns the sequence of classes in W.}
    arr:=Array(W);
    return &cat[ arr[T] : T in Keys(arr)];
end intrinsic;

/// Returns the sequence of representatives of the classes in $W$.
intrinsic Representatives(W::AlgEtQWECM)->SeqEnum[AlgEtQIdl]
{Returns the sequence of representatives of the classes in W.}
    w:=Map(W);
    return [w(c):c in Classes(W)];
end intrinsic;

/// Returns the neutral element of $W$.
intrinsic One(W::AlgEtQWECM)->AlgEtQWECMElt
{Returns the neutral element of W.}
    if not assigned W`One then
        w:=Map(W);
        W`One:=OneIdeal(Order(W))@@w;
    end if;
    return W`One;
end intrinsic;

/// Returns a random element of $W$.
intrinsic Random(W::AlgEtQWECM)->AlgEtQWECMElt
{Returns a random element of W.}
    return Random(Classes(W));
end intrinsic;

/// Returns the sequence of idempotent classes of W.
intrinsic Idempotents(W::AlgEtQWECM)->SeqEnum[AlgEtQWECMElt]
{Returns the sequence of idempotent classes of W.}
    return [x:x in Classes(W)|IsIdempotent(x)];
end intrinsic;

////////////////////////////////
// AlgEtQWECM : localizations //
////////////////////////////////

/// Given the weak equivalence class monoid of an order $R$ in the étale algebra $K$ and a prime $P$ of $R$, returns the weak equivalence class monoid of the unique overorder of $R$ which is locally equal to $R$ at $P$ and locally maximal at every other prime. This order is $R+P^k O$, where $O$ is the maximal order of $K$ and $k$ is a non-negative big-enough integer.
intrinsic Localization(W::AlgEtQWECM,P::AlgEtQIdl)->AlgEtQWECM
{Given the weak equivalence class monoid of an order R in the etale algebra K and a prime P of R, returns the weak equivalance class monoid of the unique overorder of R which is locally equal to R at P and locally maximal at every other prime. This order is R+P^kO, where O is the maximal order of K and k is a non-negative big-enough integer.}
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

/// Given a weak equivalence class monoid $W$ and a sequence of classes in $W$, returns whether the sequence is a generating set of $W$.
intrinsic IsGeneratingSet(W::AlgEtQWECM,seq::SeqEnum[AlgEtQWECMElt])->BoolElt
{Given a weak equivalence class monoid W and a sequence of classes in W, returns whether the sequence is a generating set of W.}
    cl:=Seqset(seq);
    Include(~cl,One(W));
    require forall{w:w in seq|w in W} : "The sequence does not consists of elements of W.";
    submonoid_old:=cl;
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

/// Given the weak equivalence class monoid $W$ of $R$ and a sequence of fractional $R$-ideals, returns whether the sequence represents a generating set of $W$.
intrinsic IsGeneratingSet(W::AlgEtQWECM,seq::SeqEnum[AlgEtQIdl])->BoolElt
{Given the weak equivalence class monoid W of R and a sequence of fractional R-ideals, returns whether the sequence represents a generating set of W.}
    w:=Map(W);
    return IsGeneratingSet(W,[I@@w:I in Seqset(seq)]);
end intrinsic;

/* TESTS

    printf "### Testing WkAbstract:";
	//AttachSpec("~/AlgEt/spec");
	//SetAssertions(2);
    SetClassGroupBounds("GRH");
	_<x>:=PolynomialRing(Integers());
    f:=x^4-1000*x^3-1000*x^2-1000*x-1000;
    K:=EtaleAlgebra(f);
    E:=EquationOrder(K);
    W,w:=WeakEquivalenceClassMonoidAbstract(E);
    x:=Random(W);
    xx:=x;
    SetRepresentative(x,RandomUnit(K)*Ideal(x));
    assert x eq xx;
    wk:=WKICM(E);
    wk1:=[w*RandomUnit(K):w in wk];
    assert #W eq 25;
    assert Seqset(Classes(W)) eq {x@@w : x in wk1};
    assert One(W) eq W!OneIdeal(E);
    assert IsOne(W!(RandomUnit(K)*OneIdeal(E)));
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
    assert IsOne(W!(RandomUnit(K)*OneIdeal(E)));
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
    assert IsOne(W!(RandomUnit(K)*OneIdeal(E)));
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

    L:=RationalsAsNumberField();
    K:=EtaleAlgebra([L,L,L]);
    a:=PrimitiveElement(K);
    E:=Order([a]);
    W:=WeakEquivalenceClassMonoidAbstract(E);
    O:=MaximalOrder(K);
    assert #WeakEquivalenceClassMonoidAbstract(O) eq 1;
    assert not IsMaximal(E);
    assert One(W) eq W!OneIdeal(E);
    assert IsOne(W!(RandomUnit(K)*OneIdeal(E)));
    assert IsIdempotent(W!MaximalOrder(K));
    assert #WeakEquivalenceClassMonoidAbstract(MaximalOrder(K)) eq 1;
    assert Random(W) in W;
    assert Seqset(Idempotents(W)) eq {W!T:T in OverOrders(E)};
    assert #W eq &*[#Localization(W,P) : P in SingularPrimes(E)];
    assert WeakEquivalenceClassMonoidAbstract(MaximalOrder(K)) eq Localization(W,PrimesAbove(11*E)[1]); // 11 is regular
    _:=Random(W)*Random(W);
    assert assigned W`MultiplicationTable;
    assert IsBass(E);
    assert IsGeneratingSet(W,Idempotents(W));
    _:=MultiplicationTable(W);

    SetAssertions(1);
    printf " all good!"; 
*/