function [ Miu, Sigma, P, U, S, V ] = MFGMulMF( Miu1, Sigma1, P1, U1, S1, V1, FM )
% Multiply a MFG density and a matrix Fisher density, and match the
% resulting density to a new MFG.
% See W Wang, T Lee, https://arxiv.org/abs/2003.02180, 2020
% Inputs: Miu1, Sigma1, P1, U1, S1, V1 - parameters for the MFG
%         FM - parameter for the matrix Fisher density
% Outputs: Miu, Sigma, P, U, S, V - parameters for the resulting MFG 

%% Matrix Fisher part
% parameters
[U,S,V] = psvd(U1*S1*V1'+FM);

% moments
EQ = pdf_MF_moment(diag(S));
EQQ = pdf_MF_moment2(diag(S));

%% fR1
% first moment
EQ1 = U1'*U*diag(EQ)*V'*V1;
EfR1 = vee(EQ1*S1-S1*EQ1');

% second moment
dU = U1'*U;
dV = V1'*V;
dS = dU'*S1*dV;
EfRfR1 = dU*EfRfRX(EQQ,dS)*dU';

%% fR2
EfR2 = [0;0;0];
EfRfR2 = EfRfRX(EQQ,S);

%% fR1 fR2
EfRfR12 = dU*EfRfRXS(EQQ,dS,diag(S));

%% estimation
SigmaMInv1 = diag([S1(2,2)+S1(3,3),S1(1,1)+S1(3,3),S1(1,1)+S1(2,2)]);
Sigmac1 = Sigma1-P1*SigmaMInv1*P1';
SigmaMInv = diag([S(2,2)+S(3,3),S(1,1)+S(3,3),S(1,1)+S(2,2)]);

Ex = Miu1+P1*EfR1;
Exx = Miu1*Miu1'+Miu1*EfR1'*P1'+P1*EfR1*Miu1'+P1*EfRfR1*P1'+Sigmac1;
ExfR2 = Miu1*EfR2'+P1*EfRfR12;

covxx = Exx-Ex*Ex';
covxfR2 = ExfR2-Ex*EfR2';
covfRfR2 = EfRfR2-EfR2*EfR2';

P = covxfR2*covfRfR2^-1;
Miu = Ex-P*EfR2;
Sigma = covxx-P*covxfR2'+P*SigmaMInv*P';

end


function [ E ] = EfRfRX( EQQ, X )

E(1,1) = X(3,1)^2*EQQ(4,4)+X(2,1)^2*EQQ(7,7)+X(3,2)^2*EQQ(5,5)...
    +X(2,2)^2*EQQ(8,8)+X(3,3)^2*EQQ(6,6)+X(2,3)^2*EQQ(9,9)...
    -2*X(3,2)*X(2,3)*EQQ(5,9)-2*X(2,2)*X(3,3)*EQQ(6,8);
E(2,2) = X(3,1)^2*EQQ(1,1)+X(1,1)^2*EQQ(7,7)+X(3,2)^2*EQQ(2,2)...
    +X(1,2)^2*EQQ(8,8)+X(3,3)^2*EQQ(3,3)+X(1,3)^2*EQQ(9,9)...
    -2*X(1,3)*X(3,1)*EQQ(1,9)-2*X(1,1)*X(3,3)*EQQ(3,7);
E(3,3) = X(2,1)^2*EQQ(1,1)+X(1,1)^2*EQQ(4,4)+X(2,2)^2*EQQ(2,2)...
    +X(1,2)^2*EQQ(5,5)+X(2,3)^2*EQQ(3,3)+X(1,3)^2*EQQ(6,6)...
    -2*X(2,1)*X(1,2)*EQQ(1,5)-2*X(1,1)*X(2,2)*EQQ(2,4);
E(1,2) = -X(3,1)*X(3,2)*EQQ(1,5)+X(3,1)*X(2,3)*EQQ(1,9)...
    -X(1,1)*X(2,1)*EQQ(7,7)...
    -X(3,2)*X(3,1)*EQQ(2,4)...
    -X(1,2)*X(2,2)*EQQ(8,8)+X(1,2)*X(3,3)*EQQ(6,8)...
    +X(3,3)*X(2,1)*EQQ(3,7)...
    +X(1,3)*X(3,2)*EQQ(5,9)-X(1,3)*X(2,3)*EQQ(9,9);
E(1,3) = X(2,1)*X(3,2)*EQQ(1,5)-X(2,1)*X(2,3)*EQQ(1,9)...
    -X(1,1)*X(3,1)*EQQ(4,4)...
    +X(2,2)*X(3,1)*EQQ(2,4)...
    -X(1,2)*X(3,2)*EQQ(5,5)+X(1,2)*X(2,3)*EQQ(5,9)...
    -X(2,3)*X(2,1)*EQQ(3,7)...
    +X(1,3)*X(2,2)*EQQ(6,8)-X(1,3)*X(3,3)*EQQ(6,6);
E(2,3) = -X(2,1)*X(3,1)*EQQ(1,1)+X(2,1)*X(1,3)*EQQ(1,9)...
    +X(1,1)*X(3,2)*EQQ(2,4)...
    -X(2,2)*X(3,2)*EQQ(2,2)...
    +X(1,2)*X(3,1)*EQQ(1,5)-X(1,2)*X(1,3)*EQQ(5,9)...
    +X(2,3)*X(1,1)*EQQ(3,7)-X(2,3)*X(3,3)*EQQ(3,3)...
    -X(1,3)*X(1,2)*EQQ(6,8);
E(2,1) = E(1,2);
E(3,1) = E(1,3);
E(3,2) = E(2,3);

end


function [ E ] = EfRfRXS( EQQ, X, S )

E(1,1) = -S(3)*X(2,2)*EQQ(6,8)+S(3)*X(3,3)*EQQ(6,6)...
    +S(2)*X(2,2)*EQQ(6,6)-S(2)*X(3,3)*EQQ(6,8);
E(1,2) = S(3)*X(2,1)*EQQ(3,7)-S(1)*X(2,1)*EQQ(7,7);
E(1,3) = S(2)*X(3,1)*EQQ(2,4)-S(1)*X(3,1)*EQQ(4,4);
E(2,1) = S(3)*X(1,2)*EQQ(6,8)-S(2)*X(1,2)*EQQ(8,8);
E(2,2) = -S(3)*X(1,1)*EQQ(3,7)+S(3)*X(3,3)*EQQ(3,3)...
    +S(1)*X(1,1)*EQQ(7,7)-S(1)*X(3,3)*EQQ(3,7);
E(2,3) = -S(2)*X(3,2)*EQQ(2,2)+S(1)*X(3,2)*EQQ(2,4);
E(3,1) = -S(3)*X(1,3)*EQQ(6,6)+S(2)*X(1,3)*EQQ(6,8);
E(3,2) = -S(3)*X(2,3)*EQQ(3,3)+S(1)*X(2,3)*EQQ(3,7);
E(3,3) = -S(2)*X(1,1)*EQQ(2,4)+S(2)*X(2,2)*EQQ(2,2)...
    +S(1)*X(1,1)*EQQ(4,4)-S(1)*X(2,2)*EQQ(2,4);

end

