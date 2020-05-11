function [ v, u, theta ] = logQua( q, format, checku )
% calculate the logarithm map form unit quaternions to pure quaternions
% q must be a 4-by-n or n-by-4 matrix
% if strcmp(format,'q') (default), v returns pure quaternions in the same
% dimension as q;
% if strcmp(format,'v'), v returns 3-vectors in the same dimension as q.
% [u,theta] returns the axis-angle representation, u 3-by-n, theta is 1-by-n
% if checku==true (default), check if q are unit quaternions

if ~exist('checku','var') || isempty(checku)
    checku = true;
end
if ~exist('format','var') || isempty(format)
    format = 'q';
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

% calculate log
normQv = sqrt(sum(q(2:4,:).^2));
u = q(2:4,:)./normQv;
theta = wrapToPi(atan2(normQv,q(1,:))*2);

% identity
indi = find(q(1,:)==1);
if ~isempty(indi)
    u(:,indi) = [1;0;0];
    theta(indi) = 0;
end

% format result
if strcmp(format,'q')
    v(2:4,:) = u.*theta/2;
elseif strcmp(format,'v')
    v = u.*theta;
else
    error('format must be q or v');
end

if tran
    v = v';
end

end

