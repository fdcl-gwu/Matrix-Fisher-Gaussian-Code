function [ S ] = hat( v )
% converts 3-vector to 3-by-3 skew symmetric matrix.
% if size(v,1) = 3, S returns a 3-3-size(v,2) matrix;
% if size(v,2) = 3, S returns a 3-3-size(v,1) matrix;
% otherwise returns error

if size(v,1) == 3
    N = size(v,2);
    S(1,2,:) = -reshape(v(3,:),1,1,N);
    S(1,3,:) = reshape(v(2,:),1,1,N);
    S(2,3,:) = -reshape(v(1,:),1,1,N);
elseif size(v,2) == 3
    N = size(v,1);
    S(1,2,:) = -reshape(v(:,3),1,1,N);
    S(1,3,:) = reshape(v(:,2),1,1,N);
    S(2,3,:) = -reshape(v(:,1),1,1,N);
else
    error('size along the first or second dimension of input must be 3');
end

S(2,1,:) = -S(1,2,:);
S(3,1,:) = -S(1,3,:);
S(3,2,:) = -S(2,3,:);
S(1,1,:) = 0;
S(2,2,:) = 0;
S(3,3,:) = 0;

end

