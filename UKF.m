function [ R, x, Sigma, stepT ] = UKF( gyro, Mea, sf, parameters )
% A standard implementation of UKF for attitude and gyro bias estimation.
% see J. L. Crassidis and F. L. Markley, “Unscented filtering for spacecraft
% attitude estimation,” Journal of guidance, control, and dynamics, 2003.
% Inputs: gyro - measured angular velocity
%         RMea - measured attitude
%         sf - gyroscope sampling frequency in Hz
%         parameters - a struct containing noise parameters and initial
%             conditions, etc
% Outputs: R - estimated attitude
%          x - estimated gyro bias
%          Sigma - estimated covariance matrix
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
    if ~GaussMea
        for nv = 1:nVecRef
            meaNoise(nv) = VM2Gau(meaNoise(nv));
        end
    end
else
    if ~GaussMea
        meaNoise = MF2Gau(meaNoise);
    end
end

%% initialization
% data containers
R = zeros(3,3,N);
x = zeros(3,N);
Sigma = zeros(6,6,N);
stepT = zeros(N-1,1);

% initialize
if exist('parameters','var') && isfield(parameters,'initValue')
    R(:,:,1) = parameters.initValue.U*parameters.initValue.V';
    x(:,1) = parameters.initValue.Miu;
    if GaussMea
        Sigma(:,:,1) = [parameters.initValue.RNoise,zeros(3);
            zeros(3),parameters.initValue.xNoise];
    else
        Sigma(:,:,1) = [MF2Gau(parameters.initValue.RNoise),zeros(3);
            zeros(3),parameters.initValue.xNoise];
    end
else
    R(:,:,1) = eye(3);
    x(:,1) = [0;0;0];
    Sigma(:,:,1) = [0.2^2*eye(3),zeros(3);zeros(3),0.05^2*eye(3)];
end

%% filter iteration
for n = 2:N
    tic;
    
    % propagate
    av = (gyro(:,n-1)+gyro(:,n))/2-x(:,n-1);
    [x1,w1] = GGetSigmaPoints([zeros(3,1);x(:,n-1)],Sigma(:,:,n-1));
    
    R_sigma = zeros(3,3,13);
    for ns = 13:-1:1
        if omegaLocal
            R_sigma(:,:,ns) = R(:,:,n-1)*expRot(x1(1:3,ns))*expRot((av-x1(4:6,ns))*dt);
            x1(1:3,ns) = logRot(R_sigma(:,:,end)'*R_sigma(:,:,ns),'v');
        else
            R_sigma(:,:,ns) = expRot((av-x1(4:6,ns))*dt)*R(:,:,n-1)*expRot(x1(1:3,ns));
            x1(1:3,ns) = logRot(R_sigma(:,:,end)'*R_sigma(:,:,ns),'v');
        end
    end
    
    R(:,:,n) = R_sigma(:,:,end);
    x(:,n) = x(:,n-1);
    
    Sigma(:,:,n) = sum(permute(x1-x1(:,end),[1,3,2]).*permute(x1-x1(:,end),[3,1,2]).*permute(w1,[3,1,2]),3) + ...
        [eye(3)*randomWalk^2*dt,zeros(3);zeros(3),eye(3)*biasInstability^2*dt];
    
    % update
    if rem(n,5)==2
        if meaIsVec
            % vector measurement
            vPredict = zeros(3*nVecRef,13);
            vMeaNoise = zeros(3*nVecRef,3*nVecRef);
            for nv = 1:nVecRef
                for ns = 1:13
                    if vecRefInertial
                        vPredict(3*(nv-1)+1:3*nv,ns) = R_sigma(:,:,ns)'*vRef(3*(nv-1)+1:3*nv);
                    else
                        vPredict(3*(nv-1)+1:3*nv,ns) = R_sigma(:,:,ns)*vRef(3*(nv-1)+1:3*nv);
                    end
                end
                vMeaNoise(3*(nv-1)+1:3*nv,3*(nv-1)+1:3*nv) = eye(3)*meaNoise(nv);
            end
            
            vMiu = sum(vPredict.*w1,2);
            vSigma = sum(permute(vPredict-vMiu,[1,3,2]).*permute(vPredict-vMiu,[3,1,2]).*permute(w1,[3,1,2]),3);
            Pxv = sum(permute(x1-x1(:,end),[1,3,2]).*permute(vPredict-vMiu,[3,1,2]).*permute(w1,[3,1,2]),3);
            
            K = Pxv*(vSigma+vMeaNoise)^-1;
            dx = K*(Mea(:,n)-vPredict(:,end));
        else
            % attitude measurement
            vSigma = sum(permute(x1(1:3,:),[1,3,2]).*permute(x1(1:3,:),[3,1,2]).*permute(w1,[3,1,2]),3);
            Pxv = sum(permute(x1-x1(:,end),[1,3,2]).*permute(x1(1:3,:),[3,1,2]).*permute(w1,[3,1,2]),3);
            
            if attMeaLocal
                vMeaNoise = meaNoise;
            else
                vMeaNoise = R(:,:,n)'*meaNoise*R(:,:,n);
            end
            
            K = Pxv*(vSigma+vMeaNoise)^-1;
            dx = K*logRot(R(:,:,n)'*Mea(:,:,n),'v');
        end
        
        R(:,:,n) = R(:,:,n)*expRot(dx(1:3));
        x(:,n) = x(:,n)+dx(4:6);
        Sigma(:,:,n) = Sigma(:,:,n)-K*(vSigma+vMeaNoise)*K';
    end
    
    stepT(n-1) = toc;
end

end


%% convert distributions: Monte Carlo
function [ Sigma ] = MF2Gau( S )

N = 100000;
R = pdf_MF_sampling(S,N);

v = logRot(R,'v');
Sigma = cov(v');

end


function [ sigmaSqr ] = VM2Gau( kappa )

N = 100000;
v = pdf_VM_sampling(kappa,[0;0;1],N);

sigmaSqr = mean([var(v(1,:),1),var(v(2,:),1)]);

end


%% sampling from von Mises-Fisher distribution
function x=pdf_VM_sampling(kappa,mu,N)
% simulate von Mises-Fisher distribution on \Sph^2
% K. Mardia and P. Jupp, Directional Statistics, 2000, pp. 172

%assert(abs(norm(mu)-1)<1e-5,'keyboard');

R=[mu null(mu') ];
if det(R) < 0
    R=R(:,[1 3 2]);
end

x=zeros(3,N);
for k=1:N
   theta=inv_cdf_VMF_theta(rand,kappa);
   phi=rand*2*pi;
   
   x0=[cos(theta); sin(theta)*cos(phi); sin(theta)*sin(phi)];
   x(:,k)=R*x0;    
end

end


function theta=inv_cdf_VMF_theta(C,kappa)
theta=acos(log(exp(kappa)-2*sinh(kappa)*C)/kappa);
end

