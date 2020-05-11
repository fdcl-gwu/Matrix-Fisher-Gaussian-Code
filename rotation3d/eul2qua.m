function [ q ] = eul2qua( e )
% euler angle is represneted as 3-2-1 body sequence (roll, pitch yaw)
% e is a 3-by-n or n-by-3 matrix
% q returns a 4-by-n or n-by-4 matrix depending on the size of e

% check size
if size(e,1) == 3
    tran = false;
elseif size(e,2) == 3
    e = e';
    tran = true;
else
    error('e must be of size 3-n or n-3');
end

% calculate
q(1,:) = cos(e(1,:)/2).*cos(e(2,:)/2).*cos(e(3,:)/2)+sin(e(1,:)/2).*sin(e(2,:)/2).*sin(e(3,:)/2);
q(2,:) = sin(e(1,:)/2).*cos(e(2,:)/2).*cos(e(3,:)/2)-cos(e(1,:)/2).*sin(e(2,:)/2).*sin(e(3,:)/2);
q(3,:) = cos(e(1,:)/2).*sin(e(2,:)/2).*cos(e(3,:)/2)+sin(e(1,:)/2).*cos(e(2,:)/2).*sin(e(3,:)/2);
q(4,:) = cos(e(1,:)/2).*cos(e(2,:)/2).*sin(e(3,:)/2)-sin(e(1,:)/2).*sin(e(2,:)/2).*cos(e(3,:)/2);

% format result
if tran
    q = q';
end

end

