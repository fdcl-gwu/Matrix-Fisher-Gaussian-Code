function [ R ] = expRot( S, checkSS )
% calculate the exponetial map from so(3) to SO(3)
% if size(V)=[3,3,n], treat V as n skew symmetric matrices;
% if size(V)=[3,n,1] or size(V)=[n,3,1], treat V as n 3-vectors;
% otherwise returns error.
% if size(V)=[3,3,n] and checkSS=true (default), check if V is truly skew
% symmetric.
% R is returned as a 3-3-n matrix

% check dimensions and skew-symmetric
if size(S,1)==3 && size(S,2)==3
    if ~exist('checkSS','var') || isempty(checkSS)
        checkSS = true;
    end
    if checkSS
        if ~isempty(find(S(1,1,:),1)) || ~isempty(find(S(2,2,:),1)) ||...
            ~isempty(find(S(3,3,:),1)) ||...
            ~isempty(find(S(1,2,:)+S(2,1,:),1)) ||...
            ~isempty(find(S(1,3,:)+S(3,1,:),1)) ||...
            ~isempty(find(S(2,3,:)+S(3,2,:),1))
            error('the first two dimension of input must be skew symmetric');
        end
    end
elseif size(S,1)==3 && size(S,3)==1
    S = hat(S);
elseif size(S,2)==3 && size(S,3)==1
    S = hat(S);
else
    error('Input must be of size 3-3-n, 3-n-1, or n-3-1');
end

% calculate exponential map
N = size(S,3);
R = zeros(3,3,N);
for n = 1:N
    normV = sqrt(S(1,2,n)^2+S(1,3,n)^2+S(2,3,n)^2);
    if normV==0
        R(:,:,n) = eye(3);
    else
        R(:,:,n) = eye(3)+sin(normV)/normV*S(:,:,n)+...
            (1-cos(normV))/normV^2*S(:,:,n)^2;
    end
end

end

