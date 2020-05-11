function [ e ] = rot2eul( R, checko, tran )
% R is a 3-by-3-by-n matrix
% if tran==false, e returns a 3-by-n matrix
% if tran==true, e returns a n-by-3 matrix
% if checko==true (default), check if R is orthogonal with determinant one.

if ~exist('checko','var') || isempty(checko)
    checko = true;
end
if ~exist('tran','var') || isempty(tran)
    tran = false;
end

% check size and orthogonality
if size(R,1)~=3 || size(R,2)~=3
    error('R must be of size 3-3-n');
end

orthogonalTolerance = 1e-10;
unitnessTolerance = 1e-10;
if checko
    if ~isempty(find(abs(sqrt(sum(R(:,1,:).^2))-1)>unitnessTolerance,1)) ||...
            ~isempty(find(abs(sqrt(sum(R(:,2,:).^2))-1)>unitnessTolerance,1)) ||...
            ~isempty(find(abs(sqrt(sum(R(:,3,:).^2))-1)>unitnessTolerance,1)) ||...
            ~isempty(find(abs(sum(R(:,1,:).*R(:,2,:)))>orthogonalTolerance,1)) ||...
            ~isempty(find(abs(sum(R(:,1,:).*R(:,3,:)))>orthogonalTolerance,1)) ||...
            ~isempty(find(abs(sum(R(:,2,:).*R(:,3,:)))>orthogonalTolerance,1)) ||...
            ~isempty(find(detM3(R)<0,1))
        error('R must be orthogonal matrices');
    end
end

% calculate
q = rot2qua(R);
e = qua2eul(q);

% format result
if tran
    e = e';
end

end


function [ d ] = detM3( R )

d = R(1,1,:).*R(2,2,:).*R(3,3,:)+R(1,2,:).*R(2,3,:).*R(3,1,:)+...
    R(1,3,:).*R(2,1,:).*R(3,2,:)-R(1,3,:).*R(2,2,:).*R(3,1,:)-...
    R(1,2,:).*R(2,1,:).*R(3,3,:)-R(1,1,:).*R(2,3,:).*R(3,2,:);
d = reshape(d,1,[],1);

end

