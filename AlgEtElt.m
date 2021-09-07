/* vim: set syntax=magma :*/

freeze;

/////////////////////////////////////////////////////
// Stefano Marseglia, Utrecht University, s.marseglia@uu.nl
// http://www.staff.science.uu.nl/~marse004/
/////////////////////////////////////////////////////

declare verbose ?????, 1;

/*TODO:

*/

// THIS IS A TEMPLATE

//------------
// What do we do here?
//------------

intrinsic ???(???::???) -> ???
{???}
  return ???;
end intrinsic;


//------------
// Creation and Printig for AlgEtElt
//------------

intrinsic Print(x::AlgEtElt)
{?}
    printf "%o", ?;
end intrinsic;

//------------
// Access attributes
//------------

intrinsic Parent(x::AlgEtElt) -> AlgEt
{returns the underlying etale algebra of A}
  return x`Algebra;
end intrinsic;

//------------
// Coercion
//------------

intrinsic IsCoercible(A::AlgEt, x::.) -> BoolElt, .
{Return whether x is coercible into A and the result if so}
    if Parent(x) cmpeq A then
        return true,x;
    elif Type(x) eq FldNumElt or Type(x) eq RngIntElt or Type(x) eq FldRatElt then
// TODO
        x1:=New(AlgEtElt);
        x1`Algebra:=A;
        x1`AlgAssElt:=AssAlgebra(A) ! x;
        return true, x1;
    else 
      return false,_;
    end if;
end intrinsic;

intrinsic '!'(A::AlgEt,x::AlgAssElt) -> AlgEtElt
{ given an etale algebra A and an element x of the underlying Associtive Algebra coerces x in A}
    bool,x:=IsCoercible(A,x);
    error if not bool, "the element cannot be coerced in the algebra";
    return x;
end intrinsic;

//CONTINUE FROM HERE



intrinsic Algebra(x::AlgEtElt) -> AlgEt
{returns the underlying etale algebra of A}
  return x`Algebra;
end intrinsic;

intrinsic '.'(A::AlgEt,i::RngIntElt) -> AlgEtElt
{A.i returns the ith generator of the Algebra A}
   alg:=AssociativeAlgebra(A);
   y:=New(AlgEtElt);
   y`AlgAssElt:=alg.i;
   y`Algebra:=A;
   return y;
end intrinsic;

//------------
// Operations
//------------

intrinsic '+'(x1::AlgEtElt,x2::AlgEtElt) -> AlgEtElt
{x1+x2}
   require Parent(x1) cmpeq Parent(x2): "the elements must be defined over the same algebra";
   x3:=New(AlgEtElt);
   x3`AlgAssElt:=(x1`AlgAssElt)+(x2`AlgAssElt);
   x3`Algebra:=Parent(x1);
   return x3;
end intrinsic;

intrinsic '+'(x1::.,x2::AlgEtElt) -> AlgEtElt
{x1+x2}
   bool,x1:=IsCoercible(Algebra(x2),x1);
   if bool then
        return x1 + x2;
   else
        error "x1 not coercible";
   end if;
end intrinsic;

intrinsic '+'(x1::AlgEtElt,x2::.) -> AlgEtElt
{x1+x2}
   bool,x2:=IsCoercible(Algebra(x1),x2);
   if bool then
        return x1 + x2;
   else
        error "x2 not coercible";
   end if;
end intrinsic;

intrinsic '-'(x::AlgEtElt) -> AlgEtElt
{-x}
   y:=New(AlgEtElt);
   y`AlgAssElt:=-(x`AlgAssElt);
   y`Algebra:=Parent(x);
   return y;
end intrinsic;

intrinsic '-'(x1::AlgEtElt,x2::AlgEtElt) -> AlgEtElt
{x1-x2}
   require Parent(x1) cmpeq Parent(x2): "the elements must be defined over the same algebra";
   x3:=New(AlgEtElt);
   x3`AlgAssElt:=(x1`AlgAssElt)-(x2`AlgAssElt);
   x3`Algebra:=Parent(x1);
   return x3;
end intrinsic;

intrinsic '-'(x1::.,x2::AlgEtElt) -> AlgEtElt
{x1-x2}
   bool,x1:=IsCoercible(Algebra(x2),x1);
   if bool then
        return x1 - x2;
   else
        error "x1 not coercible";
   end if;
end intrinsic;

intrinsic '-'(x1::AlgEtElt,x2::.) -> AlgEtElt
{x1-x2}
   bool,x2:=IsCoercible(Algebra(x1),x2);
   if bool then
        return x1 - x2;
   else
        error "x2 not coercible";
   end if;
end intrinsic;

intrinsic '*'(x1::AlgEtElt,x2::AlgEtElt) -> AlgEtElt
{x1*x2}
   require Parent(x1) cmpeq Parent(x2): "the elements must be defined over the same algebra";
   x3:=New(AlgEtElt);
   x3`AlgAssElt:=(x1`AlgAssElt)*(x2`AlgAssElt);
   x3`Algebra:=Parent(x1);
   return x3;
end intrinsic;

intrinsic '*'(x1::.,x2::AlgEtElt) -> AlgEtElt
{x1*x2}
   bool,x1:=IsCoercible(Algebra(x2),x1);
   if bool then
        return x1 * x2;
   else
        error "x1 not coercible";
   end if;
end intrinsic;

intrinsic '*'(x1::AlgEtElt,x2::.) -> AlgEtElt
{x1*x2}
   bool,x2:=IsCoercible(Algebra(x1),x2);
   if bool then
        return x1 * x2;
   else
        error "x2 not coercible";
   end if;
end intrinsic;

intrinsic '/'(x1::AlgEtElt,x2::AlgEtElt) -> AlgEtElt
{x1/x2}
   require Parent(x1) cmpeq Parent(x2): "the elements must be defined over the same algebra";
   x3:=New(AlgEtElt);
   x3`AlgAssElt:=(x1`AlgAssElt)*(x2`AlgAssElt)^-1;
   x3`Algebra:=Parent(x1);
   return x3;
end intrinsic;

intrinsic '/'(x1::.,x2::AlgEtElt) -> AlgEtElt
{x1/x2}
   bool,x1:=IsCoercible(Algebra(x2),x1);
   if bool then
        return x1 / x2;
   else
        error "x1 not coercible";
   end if;
end intrinsic;

intrinsic '/'(x1::AlgEtElt,x2::.) -> AlgEtElt
{x1/x2}
   bool,x2:=IsCoercible(Algebra(x1),x2);
   if bool then
        return x1 / x2;
   else
        error "x2 not coercible";
   end if;
end intrinsic;

intrinsic '^'(x1::AlgEtElt,n::RngIntElt) -> AlgEtElt
{x1^n}   
   x3:=New(AlgEtElt);
   x3`AlgAssElt:=(x1`AlgAssElt)^n;
   x3`Algebra:=Parent(x1);
   return x3;
end intrinsic;

intrinsic One(A::AlgEt) -> AlgEtElt
{the unit of A}   
   x3:=New(AlgEtElt);
   x3`AlgAssElt:=One(AssociativeAlgebra(A));
   x3`Algebra:=A;
   return x3;
end intrinsic;

intrinsic Zero(A::AlgEt) -> AlgEtElt
{the unit of A}   
   x3:=New(AlgEtElt);
   x3`AlgAssElt:=0*A.1;
   x3`Algebra:=A;
   return x3;
end intrinsic;

intrinsic 'eq'(x1::AlgEtElt,x2::AlgEtElt) -> BoolElt
{x1 eq x2}
   require Parent(x1) cmpeq Parent(x2): "the elements must be defined over the same algebra";
   return (x1`AlgAssElt) eq (x2`AlgAssElt);
end intrinsic;

intrinsic 'eq'(x1::.,x2::AlgEtElt) -> BoolElt
{x1 eq x2}
   bool,x1:=IsCoercible(Algebra(x2),x1);
   if bool then
        return x1 eq x2;
   else
        error "x1 not coercible";
   end if;
end intrinsic;

intrinsic 'eq'(x1::AlgEtElt,x2::.) -> BoolElt
{x1 eq x2}
   bool,x2:=IsCoercible(Algebra(x1),x2);
   if bool then
        return x1 eq x2;
   else
        error "x2 not coercible";
   end if;
end intrinsic;

intrinsic 'eq'(A1::AlgEt,A2::AlgEt) -> BoolElt
{A1 eq A2}
   return A1 cmpeq A2;
end intrinsic;

intrinsic Eltseq(x1::AlgEtElt) -> SeqEnum
{rational coordinates of x1}    
   return Eltseq(x1`AlgAssElt);
end intrinsic;

/* TEST

*/
