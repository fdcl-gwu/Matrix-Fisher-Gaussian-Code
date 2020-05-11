function [ R, MFG, stepT ] = MFGUnscented( gyro, RMea, parameters )
% Bayesian estimation for attitude and gyro bias using matrix
% Fisher-Gaussian distribution, with unscented moments propagation
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
Miu = parameters.xInit;
Sigma = parameters.initXNoise;
P = zeros(3);
U = parameters.RInit;
V = eye(3);
if parameters.GaussMea
    S = Gau2MF(parameters.initRNoise);
else
    S = parameters.initRNoise;
end

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
    
    % unscented transform for last step
    [xl,Rl,wl] = MFGGetSigmaPoints(Miu,Sigma,P,U,S,V);
    [xav,wav] = GGetSigmaPoints([0;0;0],eye(3)*randomWalk^2/dt);
    
    % propagate sigma points
    xp = zeros(3,13*7);
    Rp = zeros(3,3,13*7);
    wp = zeros(1,13*7);
    for i = 1:13
        for j = 1:7
            ind = 7*(i-1)+j;
            Rp(:,:,ind) = Rl(:,:,i)*expRot(((gyro(:,n-1)+gyro(:,n))...
                /2-xl(:,i)-xav(:,j))*dt);
            xp(:,ind) = xl(:,i);
            wp(ind) = wl(i)*wav(j);
        end
    end
    
    % recover prior distribution
    [Miu,Sigma,P,U,S,V] = MFGMLEAppro(xp,Rp,wp);
    Sigma = Sigma+eye(3)*biasInstability^2*dt;
    
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

