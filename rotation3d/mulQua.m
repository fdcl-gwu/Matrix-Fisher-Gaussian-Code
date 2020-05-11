function [ q ] = mulQua( q1, q2, checku )
% calculate the multiplication of two unit quaternions.
% q1 and q2 must be two 4-by-n or two n-by-4 matrices.
% q returns quaternions in the same dimensions of q1 and q2.
% If checku==true (default), chech is the input quaternions are unit.

if ~exist('checku','var') || isempty(checku)
    checku = true;
end

% check size and unitness
if size(q1,1)==4 && size(q2,1)==4
    tran = false;
elseif size(q1,2)==4 && size(q2,2)==4
    q1 = q1';
    q2 = q2';
    tran = true;
else
    error('q1 and q2 must be of sizes 4-n or n-4');
end

unitnessTolerance = 1e-10;
if checku
    if ~isempty(find(abs(sqrt(sum(q1.^2))-1)>unitnessTolerance,1))
        error('q1 must be unit quaternions');
    end
    if ~isempty(find(abs(sqrt(sum(q2.^2))-1)>unitnessTolerance,1))
        error('q2 must be unit quaternions');
    end
end

% multiplicate
q(1,:) = q1(1,:).*q2(1,:)-q1(2,:).*q2(2,:)-q1(3,:).*q2(3,:)-q1(4,:).*q2(4,:);
q(2,:) = q1(1,:).*q2(2,:)+q1(2,:).*q2(1,:)+q1(3,:).*q2(4,:)-q1(4,:).*q2(3,:);
q(3,:) = q1(1,:).*q2(3,:)-q1(2,:).*q2(4,:)+q1(3,:).*q2(1,:)+q1(4,:).*q2(2,:);
q(4,:) = q1(1,:).*q2(4,:)+q1(2,:).*q2(3,:)-q1(3,:).*q2(2,:)+q1(4,:).*q2(1,:);

if tran
    q = q';
end

end

