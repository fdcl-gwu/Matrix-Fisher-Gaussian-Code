function [ R ] = eul2rot( e )
% euler angle is represneted as 3-2-1 body sequence (roll, pitch yaw)
% e is a 3-by-n or n-by-3 matrix
% R returns a 3-by-3-by-n matrix

% check size
if size(e,2) == 3
    e = e';
elseif size(e,1) ~= 3
    error('e must be of size 3-n or n-3');
end

% calculate
R(1,1,:) = cos(e(2,:)).*cos(e(3,:));
R(1,2,:) = -cos(e(1,:)).*sin(e(3,:))+sin(e(1,:)).*sin(e(2,:)).*cos(e(3,:));
R(1,3,:) = sin(e(1,:)).*sin(e(3,:))+cos(e(1,:)).*sin(e(2,:)).*cos(e(3,:));
R(2,1,:) = cos(e(2,:)).*sin(e(3,:));
R(2,2,:) = cos(e(1,:)).*cos(e(3,:))+sin(e(1,:)).*sin(e(2,:)).*sin(e(3,:));
R(2,3,:) = -sin(e(1,:)).*cos(e(3,:))+cos(e(1,:)).*sin(e(2,:)).*sin(e(3,:));
R(3,1,:) = -sin(e(2,:));
R(3,2,:) = sin(e(1,:)).*cos(e(2,:));
R(3,3,:) = cos(e(1,:)).*cos(e(2,:));

end

