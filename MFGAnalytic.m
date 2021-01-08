function [ R, MFG, stepT ] = MFGAnalytic( gyro, Mea, sf, defQS, parameters )
% Bayesian estimation for attitude and gyro bias using matrix
% Fisher-Gaussian distribution, with analytical moments propagation
% See W Wang, T Lee, https://arxiv.org/abs/2003.02180, 2020
% Inputs: gyro - measured angular velocity
%         RMea - measured attitude
%         sf - gyroscope sampling frequency in Hz
%         defQS - true if using MFGI definition; false if using MFGB
%             definition
%         parameters - a struct containing noise parameters and initial
%             conditions, etc
% Outputs: R - estimated attitude
%          MFG - parameters for posterior MFGs
%          stepT - computation time for each discretized time step

N = size(gyro,2);
dt = 1/sf;

%% settings
% angular velocity noise parameters
if exist('parameters','var') && isfield(parameters,'omegaNoise')
    randomWalk = parameters.omegaNoise.randomWalk;
    biasInstability = parameters.omegaNoise.biasInstability;
else
    randomWalk = 10*pi/180;
    biasInstability = 500/3600*pi/180;
end

% measurement noise parameters
if exist('parameters','var') && isfield(parameters,'meaNoise')
    meaNoise = parameters.meaNoise;
else
    meaNoise = 0.2^2*eye(3);
end

% other settings
if exist('parameters','var') && isfield(parameters,'setting')
    omegaLocal = parameters.setting.omegaLocal;
    GaussMea = parameters.setting.GaussMea;
    meaIsVec = parameters.setting.meaIsVec;
else
    omegaLocal = true;
    GaussMea = true;
    meaIsVec = false;
end

if meaIsVec
    if exist('parameters','var') && isfield(parameters,'setting')
        vecRefInertial = parameters.setting.vecRefInertial;
        nVecRef = parameters.setting.nVecRef;
        vRef = parameters.setting.vRef;
    else
        vecRefInertial = true;
        nVecRef = 1;
        vRef = [0;0;1];
    end
else
    if exist('parameters','var') && isfield(parameters,'setting')
        attMeaLocal = parameters.setting.attMeaLocal;
    else
        attMeaLocal = true;
    end
end

% convert noise distributions
if meaIsVec
    if GaussMea
        for nv = 1:nVecRef
            meaNoise(nv) = Gau2VM(meaNoise(nv));
        end
    end
else
    if GaussMea
        meaNoise = Gau2MF(meaNoise);
    end
end

%% initialization
% initialize
if exist('parameters','var') && isfield(parameters,'initValue')
    Miu = -parameters.initValue.Miu;
    Sigma = parameters.initValue.xNoise;
    P = zeros(3);
    U = parameters.initValue.U;
    V = parameters.initValue.V;
    if GaussMea
        S = Gau2MF(parameters.initValue.RNoise);
    else
        S = parameters.initValue.RNoise;
    end
    S(1,1) = S(1,1)+2e-5;
    S(2,2) = S(2,2)+1e-5;
else
    Miu = [0;0;0];
    Sigma = 0.05^2*eye(3);
    P = zeros(3);
    U = eye(3);
    V = eye(3);
    S = Gau2MF(0.2^2*eye(3));
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

%% filter iteration
for n = 2:N
    try
    tic;
    
    % uncertainty propagation
    omega = (gyro(:,n-1)+gyro(:,n))/2;
    [Miu,Sigma,P,U,S,V] = MFGGyroProp(omega,Miu,Sigma,P,U,S,V,randomWalk*eye(3),biasInstability*eye(3),omegaLocal,defQS,dt);
    
    % update
    if rem(n,5)==2
        FMea = zeros(3,3);
        if meaIsVec
            if vecRefInertial
                for nv = 1:nVecRef
                    FMea = FMea + meaNoise(nv)*vRef(3*(nv-1)+1:3*nv)*Mea(3*(nv-1)+1:3*nv,n)';
                end
            else
                for nv = 1:nVecRef
                    FMea = FMea + meaNoise(nv)*Mea(3*(nv-1)+1:3*nv,n)*vRef(3*(nv-1)+1:3*nv)';
                end
            end
        else
            if attMeaLocal
                FMea = Mea(:,:,n)*meaNoise';
            else
                FMea = meaNoise'*Mea(:,:,n);
            end
        end
        
        [Miu,Sigma,P,U,S,V] = MFGMulMF(Miu,Sigma,P,U,S,V,FMea,defQS);
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
    catch
        break;
    end
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


function [ kappa ] = Gau2VM( sigmaSqr )

N = 100000;
v = randn(3,N)*sqrt(sigmaSqr)+[0;0;1];
v = v./sqrt(sum(v.^2));

rho = sqrt(sum(mean(v,2).^2));

options = optimoptions('fsolve','Algorithm','levenberg-marquardt',...
    'FunctionTolerance',1e-15,'Display','off');
kappa = fsolve(@(k) coth(k)-1/k-rho,1,options);

end

