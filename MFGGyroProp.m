function [ Miutdt, Sigmatdt, Ptdt, Utdt, Stdt, Vtdt ] = MFGGyroProp...
    ( omega, Miu, Sigma, P, U, S, V, Hu, Hv, omegaLocal, defQS, dt )

% calculate moments for the propagated density through a typical gyroscope
% kinematics model.
% See W Wang, T Lee, https://arxiv.org/abs/2003.02180, 2020
% inputs: Hu - angle random walk noise strength
%         Hv - bias random walk noise strength
%         omegaLocal - true if the angular velocity is measured in the
%             body-fixed frame; false if measued in the inertial frame
%         defQS - true if using MFGI definition; false if using MFGB
%             definition
%         dt - sampling period

s = diag(S);
[EQ,EQQ,EQQQ] = pdf_MF_moment123(s);

% E[R(t+dt)]
ERtdt = getERtdt(omega,Miu,P,U,V,Hu,EQ,dt,omegaLocal,defQS);

[Utdt,D,Vtdt] = psvd(ERtdt);
Stdt = diag(pdf_MF_M2S(diag(D),s));

% E[x(t+dt)vR(t+dt)]
ExvRtdt = getExvRtdt(omega,Miu,Sigma,P,U,s,V,Hu,EQ,EQQQ,Utdt,Stdt,Vtdt,dt,omegaLocal,defQS);

% E[v'R(t+dt)v'R(t+dt)]
EvRvRtdt = getEvRvRtdt(omega,Miu,P,U,s,V,Hu,EQQ,EQQQ,Utdt,Stdt,Vtdt,dt,omegaLocal,defQS);

%% maximum likelihood estimation
% other necessary moments
EvRvRt(1,1) = (s(2)^2+s(3)^2)*EQQ(6,6)-2*s(2)*s(3)*EQQ(6,8);
EvRvRt(2,2) = (s(1)^2+s(3)^2)*EQQ(3,3)-2*s(1)*s(3)*EQQ(3,7);
EvRvRt(3,3) = (s(1)^2+s(2)^2)*EQQ(2,2)-2*s(1)*s(2)*EQQ(2,4);

SigmaMInv = diag([s(2)+s(3),s(1)+s(3),s(1)+s(2)]);
Sigmac = Sigma-P*SigmaMInv*P';

Exxt = Sigmac+Miu*Miu'+P*EvRvRt*P';
covxxtdt = Exxt-Miu*Miu';
covxvRptdt = ExvRtdt;
covvRpvRptdt = EvRvRtdt;

% estimation
Ptdt = covxvRptdt*covvRpvRptdt^-1;
Miutdt = Miu;
Sigmatdt = covxxtdt-Ptdt*covxvRptdt'+Ptdt*(trace(Stdt)*eye(3)-Stdt)*Ptdt' + dt*(Hv*Hv)';

end

