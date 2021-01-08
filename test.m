clear;
close all;

addpath('Generate-Path');
addpath('Matrix-Fisher-Distribution');
addpath('rotation3d');

rng(1);

t = 60;
sf = 150;

% settings
parameters.setting.omegaLocal = true; % angular velocity is measured in the body-fixed frame
parameters.setting.gyroFail = false; % gyroscope does not have large unmodeled disturbance
parameters.setting.GaussMea = true; % measurement follows Gaussian distribution
parameters.setting.meaIsVec = true; % measurement are vectors
parameters.setting.vecRefInertial = true; % reference vector is fixed in the inertial frame
parameters.setting.attMeaLocal = false; % only used if meaIsVec==false
parameters.setting.nVecRef = 2; % 2 reference vecotrs
parameters.setting.vRef = [0;1;0;1;0;0]; % coordinates for reference vectors
parameters.meaNoise = [0.01,1]; % variance for vector measurements
parameters.initValue.RNoise = diag([1e10,1e10,1e10]); % initial attitude variance
parameters.initValue.xNoise = 0.1^2*eye(3); % initial bias variance
parameters.initValue.U = expRot([pi,0,0]); % initial attitude
parameters.initValue.V = eye(3);
parameters.initValue.Miu = [0.2;0.2;0.2]; % initial bias

% generate the reference motion
[gyro,Mea,RTrue,xTrue] = genTrig(t,sf,parameters);

% MEKF
[RMEKF,xMEKF,SigmaMEKF,TMEKF] = MEKF(gyro,Mea,sf,parameters);

% UKF
[RUKF,xUKF,SigmaUKF,TUKF] = UKF(gyro,Mea,sf,parameters);

% MFGI filter with analytical propagation
[RMFGIA,MFGIA,TMFGIA] = MFGAnalytic(gyro,Mea,sf,true,parameters);

% MFGI filter with unscented propagation
[RMFGIU,MFGIU,TMFGIU] = MFGUnscented(gyro,Mea,sf,true,parameters);

% MFGB filter with analytical propagation
[RMFGBA,MFGBA,TMFGBA] = MFGAnalytic(gyro,Mea,sf,false,parameters);

% MFGB filter with unscented propagation
[RMFGBU,MFGBU,TMFGBU] = MFGUnscented(gyro,Mea,sf,false,parameters);

% estimation error
error_attitude_MEKF = sqrt(sum(logRot(mulRot(invRot(RTrue),RMEKF),'v').^2));
error_attitude_UKF = sqrt(sum(logRot(mulRot(invRot(RTrue),RUKF),'v').^2));
error_attitude_MFGIA = sqrt(sum(logRot(mulRot(invRot(RTrue),RMFGIA),'v').^2));
error_attitude_MFGIU = sqrt(sum(logRot(mulRot(invRot(RTrue),RMFGIU),'v').^2));
error_attitude_MFGBA = sqrt(sum(logRot(mulRot(invRot(RTrue),RMFGBA),'v').^2));
error_attitude_MFGBU = sqrt(sum(logRot(mulRot(invRot(RTrue),RMFGBU),'v').^2));

error_bias_MEKF = sqrt(sum((xTrue-xMEKF).^2));
error_bias_UKF = sqrt(sum((xTrue-xUKF).^2));
error_bias_MFGIA = sqrt(sum((xTrue+MFGIA.Miu).^2));
error_bias_MFGIU = sqrt(sum((xTrue-MFGIU.Miu).^2));
error_bias_MFGBA = sqrt(sum((xTrue+MFGBA.Miu).^2));
error_bias_MFGBU = sqrt(sum((xTrue-MFGBU.Miu).^2));

% plot results
figure; hold on;
time = linspace(0,t,t*sf+1);
plot(time,error_attitude_MEKF*180/pi);
plot(time,error_attitude_UKF*180/pi);
plot(time,error_attitude_MFGIA*180/pi);
plot(time,error_attitude_MFGIU*180/pi);
plot(time,error_attitude_MFGBA*180/pi);
plot(time,error_attitude_MFGBU*180/pi);
xlabel('time (s)'); ylabel('attitude error (deg)');
title('Attitude Error');
legend('MEKF', 'UKF', 'MFGIA', 'MFGIU', 'MFGBA', 'MFGBU');

figure; hold on;
plot(time,error_bias_MEKF*180/pi);
plot(time,error_bias_UKF*180/pi);
plot(time,error_bias_MFGIA*180/pi);
plot(time,error_bias_MFGIU*180/pi);
plot(time,error_bias_MFGBA*180/pi);
plot(time,error_bias_MFGBU*180/pi);
xlabel('time (s)'); ylabel('bias error (deg/s)');
title('Gyro-bias Error');
legend('MEKF', 'UKF', 'MFGIA', 'MFGIU', 'MFGBA', 'MFGBU');

rmpath('Generate-Path');
rmpath('Matrix-Fisher-Distribution');
rmpath('rotation3d');


