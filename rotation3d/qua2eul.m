function [ e ] = qua2eul( q, checku )
% q is a 4-by-n or n-by-4 unit quaternion array
% e returns a 3-by-n or n-by-3 matrix depending on the size of q
% if checku==true (default), check if the quaternion is unit

if ~exist('checku','var') || isempty(checku)
    checku = true;
end

% check size and unitness
if size(q,1)==4
    tran = false;
elseif size(q,2)==4
    q = q';
    tran = true;
else
    error('q must be of size 4-n or n-4');
end

unitnessTolerance = 1e-10;
if checku
    if ~isempty(find(abs(sqrt(sum(q.^2))-1)>unitnessTolerance,1))
        error('q must be unit quaternions');
    end
end

% calculate
e(1,:) = atan2(2*(q(1,:).*q(2,:)+q(3,:).*q(4,:)), 1-2*(q(2,:).^2+q(3,:).^2));
e(2,:) = asin(2*(q(1,:).*q(3,:)-q(2,:).*q(4,:)));
e(3,:) = atan2(2*(q(1,:).*q(4,:)+q(2,:).*q(3,:)), 1-2*(q(3,:).^2+q(4,:).^2));

% gimbal lock
tolerance = 1e-10;
indn = find(q(1,:).*q(3,:)-q(2,:).*q(4,:)>0.5-tolerance);
e(1,indn) = 0;
e(2,indn) = pi/2;
e(3,indn) = -2*atan2(q(2,indn),q(1,indn));

inds = find(q(1,:).*q(3,:)-q(2,:).*q(4,:)<-0.5+tolerance);
e(1,inds) = 0;
e(2,inds) = -pi/2;
e(3,inds) = 2*atan2(q(2,inds),q(1,inds));

% format result
if tran
    e = e';
end

end

