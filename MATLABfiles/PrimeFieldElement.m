classdef PrimeFieldElement<handle
    %This class represents a prime field element. For construction it
    %receives two integers (a,p), where p is prime number and a is an element in
    %GF(p). NOTE: GF(p):={0,1,...p-1}, hence any received 'a' would be translated to this set of values by appling modulu p. The class provides implementation of
    %basic arithmetic operations above GF(p).
    
    properties
        received_a %the inserted integer
        a %an element in GF(p)
        p %the prime that characterizes the finite field above which a is taken
    end
    
    methods
        function obj = PrimeFieldElement(a,p)
            obj.p=p;
            obj.received_a=a;
            obj.a=mod(a,p);
        end
        
        function result=addition(obj,other_obj)
            %sanity check:
            if (obj.p~=other_obj.p)
                error('both elements must be aobve the same field');
            end

            %calculation of sum by defining a new element whose 'a' field is the sum result:
            add_obj=PrimeFieldElement((obj.a+other_obj.a),obj.p);
            result=add_obj.a;
        end

        function result=subtraction(obj,other_obj)
            %sanity check:
            if (obj.p~=other_obj.p)
                error('both elements must be aobve the same field');
            end

            %calculation of subtraction by defining a new element whose 'a' field is the subtraction result:
            subtract_obj=PrimeFieldElement((obj.a-other_obj.a),obj.p);
            result=subtract_obj.a;
        end

        function result=multiplication(obj,other_obj)
            %sanity check:
            if (obj.p~=other_obj.p)
                error('both elements must be aobve the same field');
            end

            %calculation of multiplication by defining a new element whose 'a' field is the mutliplication result:
            multiply_obj=PrimeFieldElement((obj.a*other_obj.a),obj.p);
            result=multiply_obj.a;
        end

        function result=inverse(obj)
            %calculation of the inverse element of obj:
            if (0==obj.a)
                error('element %d does not have an inverse above GF(%d). By default the returned value is 0\n' ,  obj.received_a,  obj.p);
            else
                % perform extended Euclidean Algorithm to find inverse:
                [d,s,t]=xgcd(obj.a,obj.p);
                if (1~=d)
                    error('the gcd of any non-zero element above a prime field with the prime defining the field should be 1');
                else
                    inverse_obj=PrimeFieldElement(s,obj.p);
                    result=inverse_obj.a;
                end
            end
            
        end

        function result=division(obj,other_obj)
            %sanity check:
            if (obj.p~=other_obj.p)
                error('both elements must be aobve the same field');
            end

            %calculation of division by defining a new element whose 'a' field is the division result:
            if (0==other_obj.a)
                error('The divider is equivalent to 0 above the operation field and division in 0 is not defined');
            else
                if (0==obj.a)
                    result=0;
                else
                    inverse_obj=PrimeFieldElement(inverse(other_obj),other_obj.p);
                    result=multiplication(obj, inverse_obj);
                end
            end
            
        end
        
    end
end
