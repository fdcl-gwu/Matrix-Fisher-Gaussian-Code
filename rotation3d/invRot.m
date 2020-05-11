function [ invR ] = invRot( R, checko )
% R must be 3-by-3-by-n matrix.
% invR returns the inverse of R(:,:,i) for each i.
% if checko==true (default), chech is the input matrices are orthogonal with
% determinant one.

if ~exist('checko','var') || isempty(checko)
    checko = true;
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
invR = permute(R,[2,1,3]);

end


function [ d ] = detM3( R )

d = R(1,1,:).*R(2,2,:).*R(3,3,:)+R(1,2,:).*R(2,3,:).*R(3,1,:)+...
    R(1,3,:).*R(2,1,:).*R(3,2,:)-R(1,3,:).*R(2,2,:).*R(3,1,:)-...
    R(1,2,:).*R(2,1,:).*R(3,3,:)-R(1,1,:).*R(2,3,:).*R(3,2,:);
d = reshape(d,1,[],1);

end

