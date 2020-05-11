function [ R ] = mulRot( R1, R2, checko )
% calcualte the multiplication of two set of rotation matrics. R1 and R2
% must be two 3-by-3-by-n matrices, or one 3-by-3-by-1 matrix and one
% 3-by-3-by-n matrix.
% R returns a 3-by-3-by-n matrix
% if checko==false (default), chech is the input matrices are orthogonal with
% determinant one.

if ~exist('checko','var') || isempty(checko)
    checko = false;
end

% check size and orthogonality
if size(R1,1)~=3 || size(R1,2)~=3
    error('R1 must be of size 3-3-n');
end
if size(R2,1)~=3 || size(R2,2)~=3
    error('R2 must be of size 3-3-n');
end
if size(R1,3)~=size(R2,3) && size(R1,3)~=1 && size(R2,3)~=1
    error('R1 and R2 must have the same number of rotation matrices, or one of them has only one matrix');
end

orthogonalTolerance = 1e-10;
unitnessTolerance = 1e-10;
if checko
    if ~isempty(find(abs(sqrt(sum(R1(:,1,:).^2))-1)>unitnessTolerance,1)) ||...
            ~isempty(find(abs(sqrt(sum(R1(:,2,:).^2))-1)>unitnessTolerance,1)) ||...
            ~isempty(find(abs(sqrt(sum(R1(:,3,:).^2))-1)>unitnessTolerance,1)) ||...
            ~isempty(find(abs(sum(R1(:,1,:).*R1(:,2,:)))>orthogonalTolerance,1)) ||...
            ~isempty(find(abs(sum(R1(:,1,:).*R1(:,3,:)))>orthogonalTolerance,1)) ||...
            ~isempty(find(abs(sum(R1(:,2,:).*R1(:,3,:)))>orthogonalTolerance,1)) ||...
            ~isempty(find(detM3(R1)<0,1))
        error('R1 must be orthogonal matrices');
    end
    if ~isempty(find(abs(sqrt(sum(R2(:,1,:).^2))-1)>unitnessTolerance,1)) ||...
            ~isempty(find(abs(sqrt(sum(R2(:,2,:).^2))-1)>unitnessTolerance,1)) ||...
            ~isempty(find(abs(sqrt(sum(R2(:,3,:).^2))-1)>unitnessTolerance,1)) ||...
            ~isempty(find(abs(sum(R2(:,1,:).*R2(:,2,:)))>orthogonalTolerance,1)) ||...
            ~isempty(find(abs(sum(R2(:,1,:).*R2(:,3,:)))>orthogonalTolerance,1)) ||...
            ~isempty(find(abs(sum(R2(:,2,:).*R2(:,3,:)))>orthogonalTolerance,1)) ||...
            ~isempty(find(detM3(R2)<0,1))
        error('R2 must be orthogonal matrices');
    end
end

% calculate
R(1,1,:) = R1(1,1,:).*R2(1,1,:)+R1(1,2,:).*R2(2,1,:)+R1(1,3,:).*R2(3,1,:);
R(1,2,:) = R1(1,1,:).*R2(1,2,:)+R1(1,2,:).*R2(2,2,:)+R1(1,3,:).*R2(3,2,:);
R(1,3,:) = R1(1,1,:).*R2(1,3,:)+R1(1,2,:).*R2(2,3,:)+R1(1,3,:).*R2(3,3,:);
R(2,1,:) = R1(2,1,:).*R2(1,1,:)+R1(2,2,:).*R2(2,1,:)+R1(2,3,:).*R2(3,1,:);
R(2,2,:) = R1(2,1,:).*R2(1,2,:)+R1(2,2,:).*R2(2,2,:)+R1(2,3,:).*R2(3,2,:);
R(2,3,:) = R1(2,1,:).*R2(1,3,:)+R1(2,2,:).*R2(2,3,:)+R1(2,3,:).*R2(3,3,:);
R(3,1,:) = R1(3,1,:).*R2(1,1,:)+R1(3,2,:).*R2(2,1,:)+R1(3,3,:).*R2(3,1,:);
R(3,2,:) = R1(3,1,:).*R2(1,2,:)+R1(3,2,:).*R2(2,2,:)+R1(3,3,:).*R2(3,2,:);
R(3,3,:) = R1(3,1,:).*R2(1,3,:)+R1(3,2,:).*R2(2,3,:)+R1(3,3,:).*R2(3,3,:);

end


function [ d ] = detM3( R )

d = R(1,1,:).*R(2,2,:).*R(3,3,:)+R(1,2,:).*R(2,3,:).*R(3,1,:)+...
    R(1,3,:).*R(2,1,:).*R(3,2,:)-R(1,3,:).*R(2,2,:).*R(3,1,:)-...
    R(1,2,:).*R(2,1,:).*R(3,3,:)-R(1,1,:).*R(2,3,:).*R(3,2,:);
d = reshape(d,1,[],1);

end

