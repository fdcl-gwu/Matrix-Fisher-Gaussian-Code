function [ invq ] = invQua( q, checku )
% q must be 4-by-n or n-by-4 matrix.
% invq returns the inverse of each quaternion, in the same dimension of q.
% if checku==true (default), chech is the input quaternions are unit

if ~exist('checku','var') || isempty(checku)
    checku = true;
end

% check size and unitness
if size(q,1)==4
    tran = false;
elseif size(q,2)==4
    q = q';
    tran = true;
elseif size(q,1)~=4
    error('q must be of size 4-n or n-4');
end

unitnessTolerance = 1e-10;
if checku
    if ~isempty(find(abs(sqrt(sum(q.^2))-1)>unitnessTolerance,1))
        error('q must be unit quaternions');
    end
end

invq = [q(1,:);-q(2:4,:)];

if tran
    invq = invq';
end

end

