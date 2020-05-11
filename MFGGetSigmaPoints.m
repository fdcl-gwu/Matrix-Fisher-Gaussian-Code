function [ x, R, w ] = MFGGetSigmaPoints( Miu, Sigma, P, U, S, V, wM, w0 )
% select sigma points from a matrix Fisher-Gaussian (MFG) distribution
% See W Wang, T Lee, https://arxiv.org/abs/2003.02180, 2020
% Inputs: Miu, Sigma, P, U, S, V - parameters for the MFG
%         wM - weight for the attitude part
%         w0 - weight at identity
% outputs: x, R - sigma points for the linear variable and attitude
%          w - weights for the sigma points

% default values
n = size(Miu,1);
if ~exist('w0','var') || isempty(w0)
    w0 = 1/((n+3)*2+1);
end
if ~exist('wM','var') || isempty(wM)
    wM = (1-w0)/(n+3)*3;
end

% declare the size
x = zeros(n,(n+3)*2+1);
R = zeros(3,3,(n+3)*2+1);
w = zeros(1,(n+3)*2+1);

%% Matrix Fisher part
% calculate sigma (use scaled normalization)
s = diag(S);
[c,dc] = pdf_MF_normal(s,1,1);

cyc = {[1,2,3],[2,3,1],[3,1,2]};
wm = cell(3,1);
cost = cell(3,1);
for i = 1:3
    if s(cyc{i}(2))+s(cyc{i}(3))>=1
        cost{i} = @(sigma)sigma+(1-sigma)*(log(c)+sum(s)-s(cyc{i}(1)))/...
            (s(cyc{i}(2))+s(cyc{i}(3)));
    else
        cost{i} = @(sigma)(sigma+(1-sigma)*(log(c)+sum(s)-s(cyc{i}(1)))+1/2)*...
            (s(cyc{i}(2))+s(cyc{i}(3)))-1/2;
    end
    wm{i} = @(sigma)1/(4*(1-cost{i}(sigma)))*...
        (1/c*(dc(cyc{i}(1))-dc(cyc{i}(2))-dc(cyc{i}(3))));
end

% minimum of sigma {cos(theta_i)>-sqrt(3)/2}
options = optimoptions('fsolve','Display','off');
sigmaMin = max(cellfun(@(eqn)fsolve(@(sigma)eqn(sigma)+sqrt(3)/2,0,options),cost));

% calculate sigma
sumwm = @(sigma)2*(wm{1}(sigma)+wm{2}(sigma)+wm{3}(sigma));
sigma = fsolve(@(sigma)sumwm(sigma)-wM,0,options);
if sigma <= sigmaMin
    sigma = sigmaMin;
    wM = sumwm(sigma);
    w0 = (1-wM)/(2*n+1);
end

% calculate sigma points and weights
Q = @(R)U'*R*V;
fR = @(R)vee(Q(R)*S-S*Q(R)');
for i = 1:3
    theta = acos(cost{i}(sigma));
    e = zeros(3,1);
    e(i) = 1;
    R(:,:,2*n+2*i-1) = U*expRot(theta*e)*V';
    R(:,:,2*n+2*i) = U*expRot(-theta*e)*V';
    x(:,2*n+2*i-1) = Miu+P*fR(R(:,:,2*n+2*i-1));
    x(:,2*n+2*i) = Miu+P*fR(R(:,:,2*n+2*i));
    w([2*n+2*i-1,2*n+2*i]) = wm{i}(sigma);
end

% calculate wM again to eliminate numerical errors
wM = sum(w(2*n+1:2*n+6));

%% Gaussian part
wG = max(1-w0-wM,1e-16);  % numerical issues when wG is around zero
SigmaTildeInv = diag([s(2)+s(3),s(1)+s(3),s(1)+s(2)]);
Sigmac = Sigma-P*SigmaTildeInv*P';
for i = 1:n
    e = zeros(n,1);
    e(i) = 1;
    y = sqrt(n/wG)*e;
    x(:,2*i-1) = sqrtm(Sigmac)*y+Miu;
    x(:,2*i) = -sqrtm(Sigmac)*y+Miu;
    R(:,:,2*i-1) = U*V';
    R(:,:,2*i) = U*V';
    w([2*i-1,2*i]) = wG/n/2;
end

%% mean point
x(:,end) = Miu;
R(:,:,end) = U*V';
w(end) = w0;

end

