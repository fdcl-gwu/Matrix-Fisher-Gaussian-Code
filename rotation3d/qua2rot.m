function [ R ] = qua2rot( q, checku )
% q is a 4-by-n or n-by-4 unit quaternion array
% R returns a 3-by-3-by-n matrix
% if checku==true (default), check if the quaternion is unit

unitnessTolerance = 1e-10;
if ~exist('checku','var') || isempty(checku)
    checku = true;
end

% chech size and if unit
if size(q,2)==4
    q = q';
elseif size(q,1)~=4
    error('one of the first two dimensions of q must be 4');
end

if checku
    if ~isempty(find(abs(sqrt(sum(q.^2))-1)>unitnessTolerance,1))
        error('quaternion q must be unit quaternions');
    end
end

% conversion
q = reshape(q,4,1,[]);
R(1,1,:) = q(1,1,:).^2+q(2,1,:).^2-q(3,1,:).^2-q(4,1,:).^2;
R(1,2,:) = 2*q(2,1,:).*q(3,1,:)-2*q(1,1,:).*q(4,1,:);
R(1,3,:) = 2*q(2,1,:).*q(4,1,:)+2*q(1,1,:).*q(3,1,:);
R(2,1,:) = 2*q(2,1,:).*q(3,1,:)+2*q(1,1,:).*q(4,1,:);
R(2,2,:) = q(1,1,:).^2-q(2,1,:).^2+q(3,1,:).^2-q(4,1,:).^2;
R(2,3,:) = 2*q(3,1,:).*q(4,1,:)-2*q(1,1,:).*q(2,1,:);
R(3,1,:) = 2*q(2,1,:).*q(4,1,:)-2*q(1,1,:).*q(3,1,:);
R(3,2,:) = 2*q(3,1,:).*q(4,1,:)+2*q(1,1,:).*q(2,1,:);
R(3,3,:) = q(1,1,:).^2-q(2,1,:).^2-q(3,1,:).^2+q(4,1,:).^2;

end

