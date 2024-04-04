classdef FiniteFieldElement<handle
    %This class represents a finite field element. Assuming the finite field, l, is an 'n' dimension extension field of a prime field GF(p),
    % the constructor is given l and a vector of length n which represents an element on the field
    % %NOTES: 
    % 1. It is assumed given vector is of the form
    %    [a_0,a_1,...,a_n-1] which equivalent to the polynomial
    %    a_0+a_1X+,...,a_n-1X^(n-1)
    % 2. Even if the given coefficients {a_0,...,a_n-1} exceeds the range [0,p-1], they would be
    %    translated to this range via mod p operation
    
    
    properties
        received_a %the given vector
        a %the given vector after converting its elements to be above GF(p) and of proper length
        l %the finite field above which a is considered
        is_0 %a boolean variable indicating whether the given a is the 0 elemnet of the field
        matrix_representation % representation of the given element a as a matrix above GLn(GF(p))
    end
    
    methods
        function obj = FiniteFieldElement(a,l)
            %examining validity of input:
            if (length(a)>l.f_x_degree)
                error('element''s degree above the given field should not exceed %d', l.f_x_degree);
            end
            %feeding properties:
            obj.l=l;
            obj.received_a=a;
            obj.a=mod(a,l.p);

            %verifying whether the element is the zero element:
            if (0==sum(obj.a))%all coefficients in obj.a are non-negative, hence their sum equals 0 only if all coefficients equal 0
                obj.is_0=1;
            else
                obj.is_0=0;
            end

            %calculating matrix representation of using the base
            %{1,x,...,x^(n-1)), where n=l.f_x_degree:
            if (1==obj.is_0)
                obj.matrix_representation=zeros(l.f_x_degree, l.f_x_degree);
            else
                obj.matrix_representation=calcMatrixRepresentationaboveFiniteField(obj.a,l);
            end
        end

        function prettyPrinting(obj)
            %printing the different representations of the element
            %as polynomial:
            fprintf('Polynomial Representation:\n');
            for i=1:obj.l.f_x_degree
                if (0~=obj.a(i))
                    if ((1<obj.a(i))||(1==i))
                        fprintf('%d',obj.a(i));
                    end
                    
                    if (2==i)
                        fprintf('x');
                    end
                    if (2<i)
                        fprintf('x^{%d}',(i-1));
                    end
                    if (0~=sum(obj.a((i+1):end)))
                        fprintf('+');
                    end
                end
                
            end
            fprintf(newline);

            %as a vector:
            fmt=['Vector Representation:' repmat(' %1.0f ',1,numel(obj.a)) '\n'];
            fprintf(fmt,obj.a);

            %as a matrix
            fprintf('Matrix Representation:\n');
            disp(obj.matrix_representation);
            newline;        
        end
        
        function result=addition(obj,other_obj)
            %sanity check:
            if (obj.l~=other_obj.l)
                error('both elements must be aobve the same field');
            end

            %calculation of sum by defining a new element whose 'a' field is the sum result:
            add_obj=FiniteFieldElement((obj.a+other_obj.a),obj.l);
            result=add_obj.a;
        end

        function result=subtraction(obj,other_obj)
            %sanity check:
            if (obj.l~=other_obj.l)
                error('both elements must be aobve the same field');
            end

            %calculation of subtraction by defining a new element whose 'a' field is the subtraction result:
            subtract_obj=FiniteFieldElement((obj.a-other_obj.a),obj.l);
            result=subtract_obj.a;
        end

        function result=multiplication(obj,other_obj)
            %sanity check:
            if (obj.l~=other_obj.l)
                error('both elements must be aobve the same field');
            end

            %calculation of multiplication by using matrix representation:
            mul_matrix=obj.matrix_representation*other_obj.matrix_representation;
            prime=obj.l.p;
            formed_element_matrix=mod(mul_matrix,prime);%matrix elements should be above the prime which serves as the kernel of the finite field
            result=(formed_element_matrix(:,1)).';
        end

        function result=division(obj,other_obj)
            %sanity check:
            if (obj.l~=other_obj.l)
                error('both elements must be aobve the same field');
            end

            %calculation of division by defining a new element whose 'a' field is the division result:
            if (1==other_obj.is_0)
                error('The divider is equivalent to 0 above the operation field and division in 0 is not defined');
            else
                if (1==obj.is_0)
                    result=zeros(1,obj.l.f_x_degree);
                else
                    %it is desired to calculate multiplication by
                    %'obj.matrix_representation*inverse(other_obj.matrix_representation)'.
                    %Only that matrix inversion may result in a fractional
                    %representation which a simple modulo wouldn't resolve.
                    %To combat this we rely on the formula:
                    %A^(-1)=det(A)^(-1)*Adjugate(A)
                    %The result of Adjugate(A) is assembled with integers (except perhaps a negligible deviation resulting from calculation is performed, hence round is used) 
                    determinant_of_divider=mod(det(other_obj.matrix_representation), other_obj.l.p);
                    inverse_determinant=inverse(PrimeFieldElement(determinant_of_divider,other_obj.l.p));
                    inverse_matrix=inverse_determinant*round(adjoint(other_obj.matrix_representation));
                    inverse_matrix=mod(inverse_matrix,other_obj.l.p);
                    division_matrix=obj.matrix_representation*inverse_matrix;
                    formed_element_matrix=mod(division_matrix,obj.l.p);
                    result=(formed_element_matrix(:,1)).';
                end
            end
            
        end
        
    end
end
