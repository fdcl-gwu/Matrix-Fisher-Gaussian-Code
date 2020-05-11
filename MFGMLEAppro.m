function [ Miu, Sigma, P, U, S, V ] = MFGMLEAppro( x, R, w, convention )
% Approximate maximum likelihood estimator for MFG. The marginal MLE for
% the attitude is solved first, which is used in solving the conditional
% MLE for the linear variable.
% See W Wang, T Lee, https://arxiv.org/abs/2003.02180, 2020
% Inputs: x, R - samples, x is N-by-Ns, R is 3-by-3-by-Ns
%         w - weights for the samples
%         convention - convention for proper singular value decomposition
% Outputs: Miu, Sigma, P, U, S, V - parameters for the estimated MFG

% default SVD convention is psvd
if ~exist('convention','var') || isempty(convention)
    convention = 1;
end

N = size(x,1);
Ns = size(R,3);

if ~exist('w','var') || isempty(w)
    w = ones(1,Ns)/Ns;
end

% empirical moments
Ex = sum(x.*w,2);
ER = sum(R.*reshape(w,1,1,[]),3);

covxx = zeros(N,N);
for ns = 1:Ns
    covxx = covxx+(x(:,ns)-Ex)*(x(:,ns)-Ex)'*w(ns);
end

% Matrix Fisher part
[U,D,V] = psvd(ER);
S = diag(pdf_MF_M2S(diag(D)));
if convention==2
    if S(3,3)<0
        S(:,1:2) = -S(:,1:2);
        U(:,1:2) = -U(:,1:2);
    end
end

% vR
Q = mulRot(mulRot(U',R,0),V,0);
vR = vee(mulRot(Q,S,0)-mulRot(S,permute(Q,[2,1,3]),0));

% other empirical moments
EvR = sum(vR.*w,2);

covxvR = zeros(N,3);
covvRvR = zeros(3,3);
for ns = 1:Ns
    covxvR = covxvR+(x(:,ns)-Ex)*(vR(:,ns)-EvR)'*w(ns);
    covvRvR = covvRvR+(vR(:,ns)-EvR)*(vR(:,ns)-EvR)'*w(ns);
end

% correlation part
P = covxvR*covvRvR^-1;

% Gaussian part
SigmaTilde2Inv = diag([S(2,2)+S(3,3),S(1,1)+S(3,3),S(1,1)+S(2,2)]);
Miu = Ex-P*EvR;
Sigma = covxx-covxvR*covvRvR^-1*covxvR'+...
    P*SigmaTilde2Inv*P';

end

