function [ R, MFG, stepT ] = MFGAnalytic( gyro, RMea, parameters )
% Bayesian estimation for attitude and gyro bias using matrix
% Fisher-Gaussian distribution, with analytical moments propagation
% See W Wang, T Lee, https://arxiv.org/abs/2003.02180, 2020
% Inputs: gyro - measured angular velocity
%         RMea - measured attitude
%         parameters - a struct containing noise parameters and initial
%             conditions, etc
% Outputs: R - estimated attitude
%          MFG - parameters for posterior MFGs
%          stepT - computation time for each discretized time step

N = size(gyro,2);
dt = parameters.dt;

% noise parameters
randomWalk = parameters.randomWalk;
biasInstability = parameters.biasInstability;
if parameters.GaussMea
    SM = Gau2MF(parameters.rotMeaNoise);
else
    SM = parameters.rotMeaNoise;
end

% initialize distribution
Miu = -parameters.xInit;
Sigma = parameters.initXNoise;
P = zeros(3);
U = parameters.RInit;
V = eye(3);
if parameters.GaussMea
    S = Gau2MF(parameters.initRNoise);
else
    S = parameters.initRNoise;
end
S(2,2) = S(2,2)+1e-5;
S(3,3) = S(3,3)+2e-5;

% data containers
MFG.Miu = zeros(3,N); MFG.Miu(:,1) = Miu;
MFG.Sigma = zeros(3,3,N); MFG.Sigma(:,:,1) = Sigma;
MFG.P = zeros(3,3,N); MFG.P(:,:,1) = P;
MFG.U = zeros(3,3,N); MFG.U(:,:,1) = U;
MFG.V = zeros(3,3,N); MFG.V(:,:,1) = V;
MFG.S = zeros(3,N); MFG.S(:,1) = diag(S);
R = zeros(3,3,N); R(:,:,1) = U*V';
stepT = zeros(N-1,1);

% filter iteration
for n = 2:N
    tic;
    % uncertainty propagation
    omega = (gyro(:,n-1)+gyro(:,n))/2;
    [Miu,Sigma,P,U,S,V] = MFGGyroProp(omega,Miu,Sigma,P,U,S,V,randomWalk*eye(3),biasInstability^2*dt*eye(3),dt);
    
    % update
    if rem(n,5)==0
        [Miu,Sigma,P,U,S,V] = MFGMulMF(Miu,Sigma,P,U,S,V,RMea(:,:,n)*SM);
    end
    
    % record results
    MFG.Miu(:,n) = Miu;
    MFG.Sigma(:,:,n) = Sigma;
    MFG.P(:,:,n) = P;
    MFG.U(:,:,n) = U;
    MFG.V(:,:,n) = V;
    MFG.S(:,n) = diag(S);
    R(:,:,n) = U*V';
    
    stepT(n-1) = toc;
end

end


function [ S ] = Gau2MF( Sigma )

N = 100000;
v = mvnrnd([0;0;0],Sigma,N);

R = expRot(v);
ER = mean(R,3);
[~,D,~] = psvd(ER);

S = diag(pdf_MF_M2S(diag(D)));

end

