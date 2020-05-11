function [ v ] = vee( S, ifTranspose, checkSS )
% convert 3-by-3 skew symmetric matrix to 3-vector
% S must be of deminsion 3-3-n, and every first two dimensions must be skew
% symmetric
% if ifTranspose=false (default), v returns a 3-n matrix
% if ifTranspose=true, v returns a n-3 matrix
% if checkSS=true (default), check if S is skew symmetric

if ~exist('ifTranspose','var') || isempty(ifTranspose)
    ifTranspose = false;
end
if ~exist('checkSS','var') || isempty(checkSS)
    checkSS = true;
end

if size(S,1)~=3 || size(S,2)~=3
    error('input must be a 3-3-n matrix');
end

tolerance = 1e-10;
if checkSS
    if ~isempty(find(abs(S(1,1,:))>tolerance,1)) ||...
            ~isempty(find(abs(S(2,2,:))>tolerance,1)) ||...
            ~isempty(find(abs(S(3,3,:))>tolerance,1)) ||...
            ~isempty(find(abs(S(1,2,:)+S(2,1,:))>tolerance,1)) ||...
            ~isempty(find(abs(S(1,3,:)+S(3,1,:))>tolerance,1)) ||...
            ~isempty(find(abs(S(2,3,:)+S(3,2,:))>tolerance,1))
        error('the first two dimension of input must be skew symmetric');
    end
end

v(1,:) = reshape(S(3,2,:),1,[]);
v(2,:) = reshape(S(1,3,:),1,[]);
v(3,:) = reshape(S(2,1,:),1,[]);

if ifTranspose
    v = v.';
end

end

