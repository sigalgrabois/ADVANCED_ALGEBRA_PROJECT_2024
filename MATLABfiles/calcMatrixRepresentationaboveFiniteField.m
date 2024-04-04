function element_matrix_representation = calcMatrixRepresentationaboveFiniteField(a, l)
%The function receives:
%'a'-an element above l, represented as a vector [a0, a1,...,a(n-1)] which
%equivalent to the polynomial a0+a1X+...a(n-1)X^(n-1)
%'l'- a finite field

%The function calculates a matrix representation of the given element, a,
%which is formed by the elementary basis {1,x,...,x^(n-1)} and the linear transformation x^i->a*x^i s.t.
%'element_matrix_representation'=[a, aX,..., aX^(n-1)]

n=l.f_x_degree;%extension dimension of the finite field
p=l.p;%the prime which used as the kernel of the extended field l
element_matrix_representation=zeros(n,n);%would hold the matrix representation of a

%initialization:
temp_polynomial=a;
element_matrix_representation(:,1)=temp_polynomial.';%first column equals the element itself   
for i=2:n
       %multiplication in x:
       %if the highest degree of 'temp_polynomial' (the last element) is
       %zero, then a cyclic shift one place right of it, is sufficient to
       %represent multiplication. Otherwise, the result should be the sum of the congruate equivalent of the highest degree and the shift right of the remaining elements of 'temp_polynomial' 
       shift_right_polynomial=[0,temp_polynomial(1:(n-1))];
       if (0~=temp_polynomial(n))%then the congruate equivalent o the highest degree should be calculated
           temp_congurant_equi=mod((temp_polynomial(n).*l.congruate_equivalency),p);%mod p is to ensure coefficients are above GF(p)
           temp_polynomial=mod((temp_congurant_equi+shift_right_polynomial),p);
       else%cyclic shift right is sufficient
           temp_polynomial=shift_right_polynomial;
       end
       element_matrix_representation(:,i)=temp_polynomial.';
end

end