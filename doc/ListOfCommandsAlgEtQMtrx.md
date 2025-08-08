<pre><b> Print</b>(x::AlgEtQMtrx)</pre>
Print the element.
<pre><b> Universe</b>(x::AlgEtQMtrx) -> AlgEtQ</pre>
Returns the étale algebra to which the entries of the matrix belongs to.
<pre><b> NumberOfRows</b>(x::AlgEtQMtrx) -> RngIntElt</pre>
Returns the number of rows of the matrix.
<pre><b> NumberOfColumns</b>(x::AlgEtQMtrx) -> RngIntElt</pre>
Returns the number of colums of the matrix.
<pre><b> Components</b>(x::AlgEtQMtrx) -> Tup</pre>
Returns the tuple of components of x. If x is defined over K=K1x...xKn with each Ki a number field, the ith compoenent is a matrix over Ki.
<pre><b> Entries</b>(x::AlgEtQMtrx) -> SeqEnum[SeqEnum[AlgEtElt]]</pre>
Returns a sequence of sequences containing the entries of the matrix.
<pre><b> Matrix</b>(entries::SeqEnum[SeqEnum[AlgEtQElt]]) -> AlgEtQMtrx</pre>
Given a sequence of s sequences of r elements of an étale algebra, it returns the corresponding sxr matrix.
<pre><b> Matrix</b>(s::RngIntElt, r::RngIntElt, Q::SeqEnum[AlgEtQElt]) -> AlgEtQMtrx</pre>
Given a sequence of sxr elements of an étale algebra, it returns the corresponding sxr matrix.
<pre><b> DiagonalMatrix</b>(Q::SeqEnum[AlgEtQElt]) -> AlgEtQMtrx</pre>
Given a sequence of elements of an étale algebra, it returns the corresponding diagonal matrix.
<pre><b> '[,]'</b>(x::AlgEtQMtrx,i::RngIntElt,j::RngIntElt) -> AlgEtQElt</pre>
Returns x[i,j].
<pre><b> Rows</b>(x::AlgEtQMtrx) -> SeqEnum[SeqEnum[AlgEtQElt]]</pre>
Returns the sequence of rows.
<pre><b> Columns</b>(x::AlgEtQMtrx) -> SeqEnum[SeqEnum[AlgEtQElt]]</pre>
Returns the sequence of columns.
<pre><b> Transpose</b>(x::AlgEtQMtrx) -> AlgEtQMtrx</pre>
Returns the transpose matrix.
<pre><b> IsSquareMatrix</b>(x::AlgEtQMtrx) -> BoolElt</pre>
Returns whether the matrix is a square matrix.
<pre><b> ScalarMatrix</b>(n::RngIntElt,s::AlgEtQElt) -> AlgEtQMtrx</pre>
Returns the scalar matrix s*I where I is the nxn identity matrix.
<pre><b> IdentityMatrix</b>(A::AlgEtQ, n::RngIntElt) -> AlgEtQMtrx</pre>
Returns the nxn identity matrix over the étale algebra A.
<pre><b> ZeroMatrix</b>(A::AlgEtQ, n::RngIntElt) -> AlgEtQMtrx</pre>
Returns the nxn identity matrix over the étale algebra A.
<pre><b> RandomMatrix</b>(A::AlgEtQ, s::RngIntElt, r::RngIntElt : bd:=3) -> AlgEtQMtrx</pre>
Returns a random sxr identity matrix over the étale algebra A, with coefficients bounded by the VarArg bd.
<pre><b> RandomMatrix</b>(A::AlgEtQ, n::RngIntElt : bd:=3) -> AlgEtQMtrx</pre>
Returns a random square nxn identity matrix over the étale algebra A, with coefficients bounded by the VarArg bd.
<pre><b> IsInvertible</b>(x::AlgEtQMtrx) -> BoolElt,AlgEtQMtrx</pre>
Returns whether the matrix is invertible and, if so, it returns also the inverse.
<pre><b> Rank</b>(x::AlgEtQMtrx) -> RngIntElt</pre>
The rank of x, which is defined as the minimum of the rank of his components.
<pre><b> Trace</b>(x::AlgEtQMtrx) -> AlgEtQElt</pre>
Determinant.
<pre><b> Determinant</b>(x::AlgEtQMtrx) -> AlgEtQElt</pre>
Determinant.
<pre><b> 'eq'</b>(x1::AlgEtQMtrx,x2::AlgEtQMtrx) -> BoolElt</pre>
Is x1=x2 ?
<pre><b> '+'</b>(x1::AlgEtQMtrx,x2::AlgEtQMtrx) -> AlgEtQMtrx</pre>
x1+x2.
<pre><b> '-'</b>(x1::AlgEtQMtrx) -> AlgEtQMtrx</pre>
-x1.
<pre><b> '-'</b>(x1::AlgEtQMtrx,x2::AlgEtQMtrx) -> AlgEtQMtrx</pre>
x1-x2.
<pre><b> '*'</b>(x1::AlgEtQMtrx,x2::AlgEtQMtrx) -> AlgEtQMtrx</pre>
x1*x2.
<pre><b> '*'</b>(x1::.,x2::AlgEtQMtrx) -> AlgEtQMtrx</pre>
x1*x2.
<pre><b> '*'</b>(x1::AlgEtQMtrx,x2::.) -> AlgEtQMtrx</pre>
x1*x2.
<pre><b> '*'</b>(x1::RngIntElt,x2::AlgEtQMtrx) -> AlgEtQMtrx</pre>
x1*x2.
<pre><b> '*'</b>(x1::FldRatElt,x2::AlgEtQMtrx) -> AlgEtQMtrx</pre>
x1*x2.
<pre><b> '*'</b>(x1::AlgEtQMtrx,x2::RngIntElt) -> AlgEtQMtrx</pre>
x1*x2.
<pre><b> '*'</b>(x1::AlgEtQMtrx,x2::FldRatElt) -> AlgEtQMtrx</pre>
x1*x2.
<pre><b> '^'</b>(x::AlgEtQMtrx,n::RngIntElt) -> AlgEtQMtrx</pre>
x^n.
<pre><b> '&+'</b>(seq::SeqEnum[AlgEtQMtrx]) -> AlgEtQMtrx</pre>
Given a sequence of AlgEtQMtrx returns the sum of the matrices.
<pre><b> '&*'</b>(seq::SeqEnum[AlgEtQMtrx]) -> AlgEtQMtrx</pre>
Given a sequence of AlgEtQMtrx returns the product of the matrices.
<pre><b> InsertBlock</b>(x::AlgEtQMtrx,y::AlgEtQMtrx,i::RngIntElt,j::RngIntElt) -> AlgEtQMtrx</pre>
Return the matrix obtained by inserting matrix Y in X at position [i,j].
<pre><b> Solution</b>(x::AlgEtQMtrx,y::AlgEtQMtrx)->AlgEtQMtrx</pre>
Returns a solution V to the system V.x=y.
<pre><b> ???</b>(x::AlgEtQMtrx) -> ???</pre>
????.
