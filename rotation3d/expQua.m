function [ q ] = expQua( v, check0 )
% calculate the exponential map from pure quaternion to unit quaternion
% If v is a 3-by-n or n-by-3 matrix, treat v as n 3-vectors
% If v is a 4-by-n or n-by-4 matrix, treat v as n pure quaternions
% in this case, if check0==true (default), check if v is pure quaternions
% q returns n unit quaternions in the same dimension as v

if ~exist('check0','var') || isempty(check0)
    check0 = true;
end

% check dimensions and pure quaternion
if size(v,1)==3
    tran = false;
    v = v/2;
elseif size(v,2)==3
    v = v'/2;
    tran = true;
elseif size(v,1)==4
    if check0
        if ~isempty(find(v(1,:)~=0,1))
            error('v must be pure quaternions');
        end
    end
    tran = false;
    v = v(2:4,:);
elseif size(v,2)==4
    if check0
        if ~isempty(find(v(:,1)~=0,1))
            error('v must be pure quaternions');
        end
    end
    tran = true;
    v = v(:,2:4)';
else
    error('v must be of size 4-n, n-4 for pure quaternions or 3-n, n-3 for vectors')
end

% calculate exponential map
theta = sqrt(sum(v.^2));
q = [cos(theta);v./theta.*sin(theta)];

if ~isempty(find(theta==0,1))
    q(:,theta==0) = [1;0;0;0];
end

% format result
if tran
    q = q';
end

end

