function [c_return,dc_return,ddc_return]=pdf_MF_normal(s,bool_scaled,bool_dc,bool_ddc)
%pdf_MF_normal: the normalizing constant for the matrix Fisher distribution
%on SO(3)
%   c = pdf_MF_normal(s) is the normalizing constant for the 
%   matrix Fisher distribution on SO(3), for a given 3x1 (or 1x3) proper singular
%   values s.
%
%   c = pdf_MF_normal(s,BOOL_SCALED) returns an exponentially scaled value 
%   specified by BOOL_SCALED:
%       0 - (defalut) is the same as pdf_MF_normal(s)
%       1 - returnes an exponentially scaled normlaizing constant,
%       exp(-sum(s))*c
%
%   See T. Lee, "Bayesian Attitude Estimation with the Matrix Fisher
%   Distribution on SO(3)", 2017, http://arxiv.org/abs/1710.03746
%
%   See also PDF_MF_NORMAL_APPROX

assert(or(min(size(s)==[1 3]),min(size(s)==[3 1])),'ERROR: s should be 3 by 1 or 1 by 3');

% if bool_scaled is not defined, then set it false
if nargin < 2
    bool_scaled=false;
end
if nargin < 3
    bool_dc=false;
end
if nargin < 4
    bool_ddc=false;
end
if bool_ddc
    bool_dc = true;
end

%% normalizing constant
if ~bool_scaled
    % return the normalizing constant without any scaling
    c = integral(@(u) f_kunze_s(u,s),-1,1);
    c_return = c;
else
    % return the normalizing constant scaled by exp(-sum(s))
    if s(1)>= s(2)
        c_bar = integral(@(u) f_kunze_s_scaled_1(u,s),-1,1);
    else
        c_bar = integral(@(u) f_kunze_s_scaled_2(u,s),-1,1);
    end
    c_return = c_bar;
end

if ~bool_dc
    return;
end

%% first order derivative
if ~bool_scaled
    dc=zeros(3,1);
    
    % derivatives of the normalizing constant
    for i=1:3
        dc(i) = integral(@(u) f_kunze_s_deriv_i(u,s,i),-1,1);
    end
    dc_return = dc;
else
    % derivatives of the scaled normalizing constant
    dc_bar = zeros(3,1);
    
    for i=1:3
        index=circshift([1 2 3],[0 4-i]);
        j=index(2);
        k=index(3);
        
        dc_bar(k) = integral(@(u) f_kunze_s_deriv_scaled(u,[s(i),s(j),s(k)]),-1,1);
    end
    dc_return = dc_bar;
end

if ~bool_ddc
    return;
end

%% second order derivative
if ~bool_scaled
    A = zeros(9,9);
    b = zeros(9,1);

    for i = 1:3
        for j = 1:3
            k = setdiff(1:3,[i,j]);
            if i==j
                if abs(s(i))~=abs(s(k(1))) && abs(s(i))~=abs(s(k(2)))
                    A(3*(i-1)+j,3*(i-1)+j) = 1;
                    b(3*(i-1)+j) = c-(-dc(i)*s(i)+dc(k(1))*s(k(1)))/(s(k(1))^2-s(i)^2)-(-dc(i)*s(i)+dc(k(2))*s(k(2)))/(s(k(2))^2-s(i)^2);
                elseif abs(s(i))~=abs(s(k(1))) && abs(s(i))==abs(s(k(2))) && s(i)~=0
                    A(3*(i-1)+j,3*(i-1)+j) = 3/2;
                    A(3*(i-1)+j,3*(i-1)+k(2)) = -1/2*sign(s(i)*s(k(2)));
                    b(3*(i-1)+j) = c-(-dc(i)*s(i)+dc(k(1))*s(k(1)))/(s(k(1))^2-s(i)^2)-dc(i)/2/s(i);
                elseif abs(s(i))~=abs(s(k(1))) && s(i)==s(k(2)) && s(i)==0
                    A(3*(i-1)+j,3*(i-1)+j) = 2;
                    b(3*(i-1)+j) = c-(-dc(i)*s(i)+dc(k(1))*s(k(1)))/(s(k(1))^2-s(i)^2);
                elseif abs(s(i))==abs(s(k(1))) && abs(s(i))~=abs(s(k(2))) && s(i)~=0
                    A(3*(i-1)+j,3*(i-1)+j) = 3/2;
                    A(3*(i-1)+j,3*(i-1)+k(1)) = -1/2*sign(s(i)*s(k(1)));
                    b(3*(i-1)+j) = c-dc(i)/2/s(i)-(-dc(i)*s(i)+dc(k(2))*s(k(2)))/(s(k(2))^2-s(i)^2);
                elseif abs(s(i))==abs(s(k(1))) && abs(s(i))==abs(s(k(2))) && s(i)~=0
                    A(3*(i-1)+j,3*(i-1)+j) = 2;
                    A(3*(i-1)+j,3*(i-1)+k(1)) = -1/2*sign(s(i)*s(k(1)));
                    A(3*(i-1)+j,3*(i-1)+k(2)) = -1/2*sign(s(i)*s(k(2)));
                    b(3*(i-1)+j) = c-dc(i)/s(i);
                elseif s(i)==s(k(1)) && abs(s(i))~=abs(s(k(2))) && s(i)==0
                    A(3*(i-1)+j,3*(i-1)+j) = 2;
                    b(3*(i-1)+j) = c-(-dc(i)*s(i)+dc(k(2))*s(k(2)))/(s(k(2))^2-s(i)^2);
                else
                    A(3*(i-1)+j,3*(i-1)+j) = 3;
                    b(3*(i-1)+j)=c;
                end
            else
                if abs(s(i))~=abs(s(j))
                    A(3*(i-1)+j,3*(i-1)+j) = 1;
                    b(3*(i-1)+j) = dc(k)+(-dc(i)*s(j)+dc(j)*s(i))/(s(j)^2-s(i)^2);
                elseif abs(s(i))==abs(s(j)) && s(i)~=0
                    A(3*(i-1)+j,3*(i-1)+j) = 3/2;
                    A(3*(i-1)+j,3*(i-1)+i) = -1/2*sign(s(i)*s(j));
                    b(3*(i-1)+j) = dc(k)-dc(j)/2/s(i);
                else
                    A(3*(i-1)+j,3*(i-1)+j) = 2;
                    b(3*(i-1)+j) = dc(k);
                end
            end
        end
    end

    ddc = A\b;
    ddc = [ddc(1:3),ddc(4:6),ddc(7:9)];
    ddc_return = ddc;
else
    A = zeros(9,9);
    b = zeros(9,1);

    for i = 1:3
        for j = 1:3
            k = setdiff(1:3,[i,j]);
            if i==j
                if abs(s(i))~=abs(s(k(1))) && abs(s(i))~=abs(s(k(2)))
                    A(3*(i-1)+j,3*(i-1)+j) = 1;
                    b(3*(i-1)+j) = -2*dc_bar(i) - c_bar/(s(i)+s(k(1)))+dc_bar(i)*s(i)/(s(k(1))^2-s(i)^2)-dc_bar(k(1))*s(k(1))/(s(k(1))^2-s(i)^2)...
                        - c_bar/(s(i)+s(k(2)))+dc_bar(i)*s(i)/(s(k(2))^2-s(i)^2)-dc_bar(k(2))*s(k(2))/(s(k(2))^2-s(i)^2);
                elseif abs(s(i))~=abs(s(k(1))) && abs(s(i))==abs(s(k(2))) && s(i)~=0
                    sig = sign(s(i)*s(k(2)));
                    A(3*(i-1)+j,3*(i-1)+j) = 3/2;
                    A(3*(i-1)+j,3*(i-1)+k(2)) = -1/2*sig;
                    b(3*(i-1)+j) = -2*dc_bar(i) - c_bar/(s(i)+s(k(1)))+dc_bar(i)*s(i)/(s(k(1))^2-s(i)^2)-dc_bar(k(1))*s(k(1))/(s(k(1))^2-s(i)^2)...
                        - (1/2-sig/2+1/2/s(i))*c_bar-(1/2/s(i)+1-sig/2)*dc_bar(i)+sig*dc_bar(k(2))/2;
                elseif abs(s(i))~=abs(s(k(1))) && s(i)==s(k(2)) && s(i)==0
                    A(3*(i-1)+j,3*(i-1)+j) = 2;
                    b(3*(i-1)+j) = -2*dc_bar(i) - c_bar/(s(i)+s(k(1)))+dc_bar(i)*s(i)/(s(k(1))^2-s(i)^2)-dc_bar(k(1))*s(k(1))/(s(k(1))^2-s(i)^2)...
                        - c_bar-2*dc_bar(i);
                elseif abs(s(i))==abs(s(k(1))) && abs(s(i))~=abs(s(k(2))) && s(i)~=0
                    sig = sign(s(i)*s(k(1)));
                    A(3*(i-1)+j,3*(i-1)+j) = 3/2;
                    A(3*(i-1)+j,3*(i-1)+k(1)) = -1/2*sign(s(i)*s(k(1)));
                    b(3*(i-1)+j) = -2*dc_bar(i) - (1/2-sig/2+1/2/s(i))*c_bar-(1/2/s(i)+1-sig/2)*dc_bar(i)+sig*dc_bar(k(1))/2 ...
                        - c_bar/(s(i)+s(k(2)))+dc_bar(i)*s(i)/(s(k(2))^2-s(i)^2)-dc_bar(k(2))*s(k(2))/(s(k(2))^2-s(i)^2);
                elseif abs(s(i))==abs(s(k(1))) && abs(s(i))==abs(s(k(2))) && s(i)~=0
                    sig1 = sign(s(i)*s(k(1)));
                    sig2 = sign(s(i)*s(k(2)));
                    A(3*(i-1)+j,3*(i-1)+j) = 2;
                    A(3*(i-1)+j,3*(i-1)+k(1)) = -1/2*sig1;
                    A(3*(i-1)+j,3*(i-1)+k(2)) = -1/2*sig2;
                    b(3*(i-1)+j) = -2*dc_bar(i) - (1/2-sig1/2+1/2/s(i))*c_bar-(1/2/s(i)+1-sig1/2)*dc_bar(i)+sig1*dc_bar(k(1))/2 ...
                        - (1/2-sig2/2+1/2/s(i))*c_bar-(1/2/s(i)+1-sig2/2)*dc_bar(i)+sig2*dc_bar(k(2))/2;
                elseif s(i)==s(k(1)) && abs(s(i))~=abs(s(k(2))) && s(i)==0
                    A(3*(i-1)+j,3*(i-1)+j) = 2;
                    b(3*(i-1)+j) = -2*dc_bar(i) - c_bar-2*dc_bar(i)...
                        - c_bar/(s(i)+s(k(2)))+dc_bar(i)*s(i)/(s(k(2))^2-s(i)^2)-dc_bar(k(2))*s(k(2))/(s(k(2))^2-s(i)^2);
                else
                    A(3*(i-1)+j,3*(i-1)+j) = 3;
                    b(3*(i-1)+j) = -2*dc_bar(i) - c_bar-2*dc_bar(i) - c_bar-2*dc_bar(i);
                end
            else
                if abs(s(i))~=abs(s(j))
                    A(3*(i-1)+j,3*(i-1)+j) = 1;
                    b(3*(i-1)+j) = -c_bar/(s(i)+s(j)) - (1+s(j)/(s(j)^2-s(i)^2))*dc_bar(i)...
                        - (1-s(i)/(s(j)^2-s(i)^2))*dc_bar(j) + dc_bar(k);
                elseif abs(s(i))==abs(s(j)) && s(i)~=0
                    sig = sign(s(i)*s(j));
                    A(3*(i-1)+j,3*(i-1)+j) = 3/2;
                    A(3*(i-1)+j,3*(i-1)+i) = -1/2*sig;
                    b(3*(i-1)+j) = -(1/2-sig/2+1/2/s(i))*c_bar - (3/2-sig)*dc_bar(i) - (3/2+1/2/s(i))*dc_bar(j) + dc_bar(k);
                else
                    A(3*(i-1)+j,3*(i-1)+j) = 2;
                    b(3*(i-1)+j) = -c_bar - 2*dc_bar(i) - 2*dc_bar(j) + dc_bar(k);
                end
            end
        end
    end

    ddc_bar = A\b;
    ddc_bar = [ddc_bar(1:3),ddc_bar(4:6),ddc_bar(7:9)];
    ddc_return = ddc_bar;
end

end


function Y=f_kunze_s(u,s)
% integrand for the normalizing constant

J=besseli(0,1/2*(s(1)-s(2))*(1-u)).*besseli(0,1/2*(s(1)+s(2))*(1+u));
Y=1/2*exp(s(3)*u).*J;

end


function Y=f_kunze_s_scaled_1(u,s)
% integrand for the normalizing constant scaled by exp(-sum(s)) when s(1)
% >= s(2)

J=besseli(0,1/2*(s(1)-s(2))*(1-u),1).*besseli(0,1/2*(s(1)+s(2))*(1+u),1);
Y=1/2*exp((s(2)+s(3))*(u-1)).*J;

end


function Y=f_kunze_s_scaled_2(u,s)
% integrand for the normalizing constant scaled by exp(-sum(s)) when s(1)
% <= s(2)

J=besseli(0,1/2*(s(1)-s(2))*(1-u),1).*besseli(0,1/2*(s(1)+s(2))*(1+u),1);
Y=1/2*exp((s(1)+s(3))*(u-1)).*J;

end


function Y=f_kunze_s_deriv_scaled(u,s)
% integrand for the derivative of the scaled normalizing constant

J=besseli(0,1/2*(s(1)-s(2))*(1-u),1).*besseli(0,1/2*(s(1)+s(2))*(1+u),1);
Y=1/2*J.*(u-1).*exp((min([s(1) s(2)])+s(3))*(u-1));

end


function Y=f_kunze_s_deriv_i(u,s,i)
% integrand for the derivative of the normalizing constant
index=circshift([1 2 3],[0 4-i]);
j=index(2);
k=index(3);

J00=besseli(0,1/2*(s(j)-s(k))*(1-u)).*besseli(0,1/2*(s(j)+s(k))*(1+u));
Y=1/2*J00.*u.*exp(s(i)*u);

end
