classdef FiniteField<handle
    %This class represents a finite field formed as an extension of a prime field, GF(p) by an irreduciable polynomial f(x). 
    % The class assumes p is indeed prime and verifies irreduciability for
    % polynomials whose degree is up to 3. It is assumed polynomial
    % coefficients are arranged from free element (leftest coefficient) to
    % the highest degree
    
    properties
        p %prime whose corresponding field is the kernel of the described finite field
        f_x_original % given irreduciable polynomial whose highest degree coefficient may be larger than 1
        f_x_monic % monic representation of the given irreduciable polynomial
        f_x_degree %the degree of the given irreduciable polynomial which also equivalent to the field extension dimension
        field_size % number of elements above the described finite field
        congruate_equivalency % the equivalence polynomial representation of x^(f_x_degree) deduced from the irreduciable polynomial
        
    end
    
    methods
        function obj = FiniteField(p,f_x)
            obj.p=p;
            obj.f_x_original=f_x;
            obj.f_x_degree=length(f_x)-1;
            %verifying irreduciability:
            if (1>=obj.f_x_degree) || (0==obj.f_x_original(1))% for irreduciability, deg(f(x)) must be at least 2 and free coefficient must be non-zero
                error('polynomial is not irreduciable, finite field can not be formed');
            else
                if (obj.f_x_degree<=3)%irreducability can be verified by root examination above the field
                    irreduce=false;
                    i=1;
                    while ((i<obj.p)&&(false==irreduce))
                        if (2==obj.f_x_degree)
                            power_vec=mod([1,i,i^2],obj.p);
                        else %i.e. deg(f(x))=3
                            power_vec=mod([1,i,i^2,i^3],obj.p);
                        end
                        test_result=mod(sum(power_vec.*obj.f_x_original),obj.p);
                        if (0==test_result)%meaning i is a root of the given polynomial hence f(x) is reducable and can't be used to build an extension field
                          irreduce=true;
                        else%moving to examine the next element above GF(p)
                            i=i+1;
                        end
                    end
                    if (true==irreduce)
                        error('polynomial is not irreduciable, finite field can not be formed');
                    end
                else %ireeducability examination is beyond this project's scope
                    warning('It is assumed the given polynomial is irreduciable. This assumption was not verified');
                end
            end

            %converting f(x) to its monic form:
            leading_coeff=obj.f_x_original(length(f_x));%the coefficient of the given polynomial's highest degree
            if (1==leading_coeff)%i.e. the given polynomial is in its monic form
                obj.f_x_monic=f_x;
            else
                %finding inverse to 'leading_coeff':
                element=PrimeFieldElement(leading_coeff,obj.p);
                inverse_coeff=inverse(element);
                obj.f_x_monic=mod((inverse_coeff.*f(x)),obj.p);
            end

            %calculating field size:
            obj.field_size=obj.p^obj.f_x_degree;%field size if p^(deg(f(x))

            %congurate equivalency: given polynomial
            %f(x)=a_0+a_1X+...+a_(n-1)X^(n-1)+X^n. perfoming modulu f(x)
            %means X^n=-a_0-a_1X-...-a_(n-1)X^(n-1)
            lower_power_coeff_vec=obj.f_x_monic(1:(length(f_x)-1));%vector of f(x)_monic coefficients without the coefficient of X^(deg(f(x))
            obj.congruate_equivalency=mod((-lower_power_coeff_vec),obj.p);
        end
        
    end
end
