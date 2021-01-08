function [ gyroMea, Mea, RTrue, biasTrue ] = genRotAxis( t, sf, parameters )
% generate a reference rotational motion with noisy angular velocity and
% attitude/vector measurements. The rotational motion is that the rigid
% body rotates about the body-fixed z-axis at 6 rad/s, and at the same time
% rotates about the inertial y-axis at 1 rad/s.
% The angular velocity has two noises, a white noise and a bias modeled as
% a Brownian motion. The attitude measurement noise is modeled as either a
% matrix Fisher distribution or a Gaussian distribution in the tangent
% space of the true attitude. The vector measurement noise is models as
% either a von Mises Fisher distribution or a Gaussian distribution.
% Input: t - total time of the trajectory
%        sf - inverse of the discretized time step
%        parameters - a struct containing noise strengthes, etc
% Output: gyroMea - measured angular velocity
%         RMea - measured attitude
%         RTrue - true attitude
%         biasTrue - true bias

time = (0:1/sf:t)';
N = length(time);

%% parameters
% motion parameters
E.initRA = [0;0;1];      % fixed
E.av = 6;
E.RARA = [0;1;0];        % fixed
E.avRA = 1;

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

if exist('parameters','var') && isfield(parameters,'setting')
    gyroFail = parameters.setting.gyroFail;
else
    gyroFail = false;
end

%% true attitude
R1 = @(t)[cos(t*E.avRA), 0, sin(t*E.avRA);
          0, 1, 0;
          -sin(t*E.avRA), 0, cos(t*E.avRA)];
R2 = @(t)[cos(t*E.av), -sin(t*E.av), 0;
          sin(t*E.av), cos(t*E.av), 0;
          0, 0, 1];
R = @(t)R1(t)*R2(t);

RTrue = zeros(3,3,N);
for n = 1:N
    RTrue(:,:,n) = R(time(n));
end
    
%% true angular velocity
dR11 = @(t)-E.avRA*sin(t*E.avRA)*cos(t*E.av)-E.av*cos(t*E.avRA)*sin(t*E.av);
dR12 = @(t)E.avRA*sin(t*E.avRA)*sin(t*E.av)-E.av*cos(t*E.avRA)*cos(t*E.av);
dR13 = @(t)E.avRA*cos(t*E.avRA);
dR21 = @(t)E.av*cos(t*E.av);
dR22 = @(t)-E.av*sin(t*E.av);
dR23 = @(t)0;
dR31 = @(t)-E.avRA*cos(t*E.avRA)*cos(t*E.av)+E.av*sin(t*E.avRA)*sin(t*E.av);
dR32 = @(t)E.avRA*cos(t*E.avRA)*sin(t*E.av)+E.av*sin(t*E.avRA)*cos(t*E.av);
dR33 = @(t)-E.avRA*sin(t*E.avRA);
dR = @(t)[dR11(t),dR12(t),dR13(t);dR21(t),dR22(t),dR23(t);dR31(t),dR32(t),dR33(t)];

gyro = zeros(3,N);
for n = 1:N
    if omegaLocal
        gyro(:,n) = vee(R(time(n))'*dR(time(n)));
    else
        gyro(:,n) = vee(dR(time(n))*R(time(n))');
    end
end

%% add noise
biasNoise = randn(3,N)*biasInstability*sqrt(sf);
biasTrue = cumsum(biasNoise/sf,2);

gyroNoise = randn(3,N)*randomWalk*sqrt(sf);
gyroMea = gyro+gyroNoise;

%% gyroscope failure
if gyroFail
    gyroMea(:,20*sf+1:22*sf) = zeros(3,2*sf);
end

%% measurement
if meaIsVec
    Mea = zeros(3*nVecRef,N);
    if GaussMea
        vecNoise = zeros(3*nVecRef,N);
        for nv = 1:nVecRef
            vecNoise(3*(nv-1)+1:3*nv,:) = randn(3,N)*sqrt(meaNoise(nv));
        end
        
        for n = 1:N
            for nv = 1:nVecRef
                if vecRefInertial
                    Mea(3*(nv-1)+1:3*nv,n) = RTrue(:,:,n)'*vRef(3*(nv-1)+1:3*nv)...
                        + vecNoise(3*(nv-1)+1:3*nv,n);
                else
                    Mea(3*(nv-1)+1:3*nv,n) = RTrue(:,:,n)*vRef(3*(nv-1)+1:3*nv)...
                        + vecNoise(3*(nv-1)+1:3*nv,n);
                end
                Mea(3*(nv-1)+1:3*nv,n) = Mea(3*(nv-1)+1:3*nv,n)./sqrt(sum(Mea(3*(nv-1)+1:3*nv,n).^2));
            end
        end
    else
        for n = 1:N
            for nv = 1:nVecRef
                if vecRefInertial
                    Mea(3*(nv-1)+1:3*nv,n) = pdf_VM_sampling(meaNoise(nv),RTrue(:,:,n)'*vRef(3*(nv-1)+1:3*nv),1);
                else
                    Mea(3*(nv-1)+1:3*nv,n) = pdf_VM_sampling(meaNoise(nv),RTrue(:,:,n)*vRef(3*(nv-1)+1:3*nv),1);
                end
            end
        end
    end
else
    if GaussMea
        RNoise = expRot(mvnrnd([0;0;0],meaNoise,N));
    else
        RNoise = pdf_MF_sampling(meaNoise,N);
    end
    
    if attMeaLocal
        Mea = mulRot(RTrue,RNoise);
    else
        Mea = mulRot(RNoise,RTrue);
    end
end

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

