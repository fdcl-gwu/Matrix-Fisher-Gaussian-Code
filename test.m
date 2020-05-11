clear;
close all;

addpath('Generate-Path');
addpath('Matrix-Fisher-Distribution');
addpath('rotation3d');

rng(2);

t = 60;
sf = 150;

% general settings
parameters = [];
parameters.t = t;
parameters.dt = 1/sf;
parameters.randomWalk = 10*pi/180;
parameters.biasInstability = 500/3600*pi/180;
parameters.GaussMea = false;
parameters.rotMeaNoise = diag([12,12,12]);
parameters.initRNoise = diag([200,200,200]);
parameters.initXNoise = 0.1^2*eye(3);
parameters.RInit = eye(3);
parameters.xInit = [0;0;0];

[gyro,RMea,RTrue,xTrue] = genTrig(t,sf,parameters);
% stachastic settings
parameters.RInit = RTrue(:,:,1)*expRot([pi,0,0]);
parameters.xInit = [0.2;0.2;0.2];

% estimation
[RMEKF,xMEKF,G,TMEKF] = MEKF(gyro,RMea,parameters);
[RMFGA,MFGA,TMFGA] = MFGAnalytic(gyro,RMea,parameters);
[RMFGU,MFGU,TMFGU] = MFGUnscented(gyro,RMea,parameters);

% estimation error
error_attitude_MFGA = sqrt(sum(logRot(mulRot(invRot(RTrue),RMFGA),'v').^2));
error_attitude_MFGU = sqrt(sum(logRot(mulRot(invRot(RTrue),RMFGU),'v').^2));
error_attitude_MEKF = sqrt(sum(logRot(mulRot(invRot(RTrue),RMEKF),'v').^2));
error_bias_MFGA = sqrt(sum((xTrue+MFGA.Miu).^2));
error_bias_MFGU = sqrt(sum((xTrue+MFGU.Miu).^2));
error_bias_MEKF = sqrt(sum((xTrue-xMEKF).^2));

% plot results
figure; hold on;
time = linspace(0,parameters.t,parameters.t*sf+1);
plot(time,error_attitude_MFGA*180/pi);
plot(time,error_attitude_MFGU*180/pi);
plot(time,error_attitude_MEKF*180/pi);
xlabel('time (s)'); ylabel('attitude error (deg)');
title('Attitude Error');
legend('MFG-Analytic','MFG-Unscented','MEKF');

figure; hold on;
plot(time,error_bias_MFGA*180/pi);
plot(time,error_bias_MFGU*180/pi);
plot(time,error_bias_MEKF*180/pi);
xlabel('time (s)'); ylabel('bias error (deg/s)');
title('Gyro-bias Error');
legend('MFG-Analytic','MFG-Unscented','MEKF');

rmpath('Generate-Path');
rmpath('Matrix-Fisher-Distribution');
rmpath('rotation3d');


