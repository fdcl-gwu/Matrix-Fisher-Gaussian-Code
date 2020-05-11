function [ x, w ] = GGetSigmaPoints( Miu, Sigma, w0 )
% select sigma points from a Gaussian vector
% Input: Miu - mean
%        Sigma - covariance matrix
%        w0 - weight at the mean
% Output: x - sigma points, a n-by-2n+1 matrix
%         w - weights for the sigma points

n = size(Miu,1);
if ~exist('w0','var') || isempty(w0)
    w0 = 1/(2*n+1);
end

wG = 1-w0;
x = zeros(3,2*n+1);
w = zeros(1,2*n+1);
for i = 1:n
    e = zeros(n,1);
    e(i) = 1;
    y = sqrt(n/wG)*e;
    x(:,2*i-1) = sqrtm(Sigma)*y+Miu;
    x(:,2*i) = -sqrtm(Sigma)*y+Miu;
    w([2*i-1,2*i]) = wG/n/2;
end

x(:,end) = Miu;
w(end) = w0;

end

