function [ att ] = int0M( omega, att0, dt, method )
% zeroth order angular velocity integration for multiple time steps
% omega must be 3-by-n or n-by-3 matrix
% att0 is either 4-by-1, 1-by-4 as a quaternion; or 3-by-3 as a rotation
% matrix. If att0 is a quaternion, att returns a 4-by-n matrix; if att0 is
% a rotation matrix, att returns a 3-by-3-by-n matrix.
% if method==0 (default), use midpoint integration;
%    method==1, use forward integration;
%    method==2, use backward integration.

if ~exist('method','var') || isempth(method)
    method = 0;
end

% check dimensions
if size(omega,2)==3
    omega = omega';
elseif size(omega,1)~=3
    error('one of the first two dimensions of omega must be of size 3');
end
if ~isscalar(dt)
    error('dt must be a scalar');
end
if size(att0,1)==4
    isQuat = true;
elseif size(att0,2)==4
    isQuat = true;
    att0 = att0';
elseif size(att0,1)==3 && size(att0,2)==3
    isQuat = false;
else
    error('att0 must be a quaternion (size 4-1,1-4) or a rotation matrix (size 3-3)')
end

% choose angular velocity for different integration methods
N = size(omega,2);
if method==0
    omega = (omega(:,1:end-1)+omega(:,2:end))/2;
elseif method==1
    omega = omega(:,1:end-1);
elseif method==2
    omega = omega(:,2:end);
else
    error('method must be 0 (midpoint), 1 (forward), or 2 (backward)');
end

% integration
if isQuat
    dQ = expQua(omega*dt);
    att = zeros(4,N);
    att(:,1) = att0;
    for n = 2:N
        att(:,n) = mulQua(att(:,n-1),dQ(:,n-1));
    end
else
    dR = expRot(omega*dt);
    att = zeros(3,3,N);
    att(:,:,1) = att0;
    for n = 2:N
        att(:,:,n) = att(:,:,n-1)*dR(:,:,n-1);
    end
end

end

