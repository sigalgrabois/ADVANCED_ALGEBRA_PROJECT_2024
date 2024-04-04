function [d,s,t]=xgcd(a,b)
%This function performs the extended Euclides algorithm which receives two
%complete numbers a and b and returns their g.c.d, (a,b)=d and two complete
%numbers s,t which uphold: d=s*a+t*b

%Sanity checks:
if ((0~=mod(a,1)) || (0~=mod(b,1)))
    error('both inputs must be complete numbers');
end

%Edge cases:
if ((0==a)&&(0==b))
    error('at least one input must not be equal to 0');
end
%The g.c.d. calculation is done for the absolute values of a,b and at the
%end result the signs are embeded at the found s,t
mag_a=abs(a);
mag_b=abs(b);
if (0==a)||(0==b)%the case of both equals to 0 was already adressed
    d=max(abs(a),abs(b));
    s=1;
    t=1;
else%an iterative process is required
    %required data structures: 
    max_stages=mod(max(mag_a,mag_b),min(mag_a,mag_b));%the first residue bounds the possible number of stages
    calc_residue_matrix=zeros(max_stages,4);%each raw represents a stage in the process where raw i represent the equation r_{i-2}=r_{i-1}*q_{i}+r_{i}. Notice that r_{-1}=a; r_{0}=b; hence cell (i,1)=r_{i-2}; cell(i,2)=r_{i-1}; cell (i,3)=q_{i}; cell (i,4)=r_{i}
    iteration_idx=0;%counts the number of stages (which at most can be min(|a|,|b|)), the procedure is finished when r_{i}=0;
    calc_residue_coefficients_matrix=zeros(max_stages,2); % raw i represent the equation: r_{iteration_idx-1}=s_{i-2}*r_{i-2}+s_{i-1}*r_{i-1}. cell (i,1)=s_{i-2};cell (i,2)=s_{i-1}
    temp_residue=1;%a default value to ensure entering the loop:
    earlier_residue=mag_a;
    early_residue=mag_b;
    while (0~=temp_residue)%a loop calculating the residues
        iteration_idx=iteration_idx+1;
        calc_residue_matrix(iteration_idx,1)=earlier_residue;
        calc_residue_matrix(iteration_idx,2)=early_residue;
        calc_residue_matrix(iteration_idx,3)=floor(earlier_residue/early_residue);
        calc_residue_matrix(iteration_idx,4)=earlier_residue-(early_residue*calc_residue_matrix(iteration_idx,3));
        %forwarding variables:
        temp_residue=calc_residue_matrix(iteration_idx,4);
        earlier_residue=early_residue;
        early_residue=temp_residue;
    end
    d=calc_residue_matrix(iteration_idx-1,4);%the last non-zero residue is the g.c.d.
    %initialize conditions:
    calc_residue_coefficients_matrix(iteration_idx,1)=0;
    calc_residue_coefficients_matrix(iteration_idx,2)=1;
    for i=(iteration_idx-1):(-1):1
        calc_residue_coefficients_matrix(i,1)=calc_residue_coefficients_matrix(i+1,2);
        calc_residue_coefficients_matrix(i,2)=calc_residue_coefficients_matrix(i+1,1)-calc_residue_coefficients_matrix(i+1,2)*calc_residue_matrix(i,3);                  
    end
    s=calc_residue_coefficients_matrix(1,1);
    t=calc_residue_coefficients_matrix(1,2);
end
%updating the coefficients signs if necessary:
s=sign(a)*s;
t=sign(b)*t;
%sanity check:
if ((d<=0)||(d~=(s*a+t*b)))
    error('something went wrong with the calculations');
end