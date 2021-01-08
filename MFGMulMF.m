function [ MiuC, SigmaC, PC, UC, SC, VC ] = MFGMulMF( Miu, Sigma, P, U, S, V, FM, defQS )
% Multiply a MFG density and a matrix Fisher density, and match the
% resulting density to a new MFG.
% See W Wang, T Lee, https://arxiv.org/abs/2003.02180, 2020
% Inputs: Miu, Sigma, P, U, S, V - parameters for the MFG
%         FM - parameter for the matrix Fisher density
%         defQS - true if using MFGI definition; false if using MFGB
%             definition
% Outputs: MiuC, SigmaC, PC, UC, SC, VC - parameters for the resulting MFG 

%% Matrix Fisher part
% parameters
[UC,SC,VC] = psvd(U*S*V'+FM);

% moments
EQC = diag(pdf_MF_moment(diag(SC)));
EQQC = pdf_MF_moment2(diag(SC));

%% vR
UT = U'*UC;
VT = V'*VC;
ST = UT'*S*VT;

% EvR
if defQS
    EvR = UT*vee(EQC*ST'-ST*EQC');
else
    EvR = VT*vee(ST'*EQC-EQC'*ST);
end

% EvRvR
if defQS
    EvRvR = UT*EvRvRST(EQQC,ST,defQS)*UT';
else
    EvRvR = VT*EvRvRST(EQQC,ST,defQS)*VT';
end

%% vRC
EvRC = [0;0;0];

EvRCvRC(1,1) = (SC(2,2)^2+SC(3,3)^2)*EQQC(6,6)-2*SC(2,2)*SC(3,3)*EQQC(6,8);
EvRCvRC(2,2) = (SC(1,1)^2+SC(3,3)^2)*EQQC(3,3)-2*SC(1,1)*SC(3,3)*EQQC(3,7);
EvRCvRC(3,3) = (SC(1,1)^2+SC(2,2)^2)*EQQC(2,2)-2*SC(1,1)*SC(2,2)*EQQC(2,4);

%% vR vRC
if defQS
    EvRvRC = UT*EfRfRXS(EQQC,ST,diag(SC),defQS);
else
    EvRvRC = VT*EfRfRXS(EQQC,ST,diag(SC),defQS);
end

%% estimation
SigmaMInv = diag([S(2,2)+S(3,3),S(1,1)+S(3,3),S(1,1)+S(2,2)]);
Sigmac = Sigma-P*SigmaMInv*P';
SigmaMInvC = diag([SC(2,2)+SC(3,3),SC(1,1)+SC(3,3),SC(1,1)+SC(2,2)]);

ExC = Miu+P*EvR;
ExxC = Miu*Miu'+Miu*EvR'*P'+P*EvR*Miu'+P*EvRvR*P'+Sigmac;
ExvRC = Miu*EvRC'+P*EvRvRC;

covxxC = ExxC-ExC*ExC';
covxvRC = ExvRC-ExC*EvRC';
covvRvRC = EvRCvRC-EvRC*EvRC';

PC = covxvRC*covvRvRC^-1;
MiuC = ExC-PC*EvRC;
SigmaC = covxxC-PC*covxvRC'+PC*SigmaMInvC*PC';

end


function [ E ] = EvRvRST( EQQ, ST, defQS )

if defQS
    E(1,1) = ST(3,1)^2*EQQ(4,4)+ST(2,1)^2*EQQ(7,7)+ST(3,2)^2*EQQ(5,5)...
        +ST(2,2)^2*EQQ(8,8)+ST(3,3)^2*EQQ(6,6)+ST(2,3)^2*EQQ(9,9)...
        -2*ST(3,2)*ST(2,3)*EQQ(5,9)-2*ST(2,2)*ST(3,3)*EQQ(6,8);
    E(2,2) = ST(3,1)^2*EQQ(1,1)+ST(1,1)^2*EQQ(7,7)+ST(3,2)^2*EQQ(2,2)...
        +ST(1,2)^2*EQQ(8,8)+ST(3,3)^2*EQQ(3,3)+ST(1,3)^2*EQQ(9,9)...
        -2*ST(1,3)*ST(3,1)*EQQ(1,9)-2*ST(1,1)*ST(3,3)*EQQ(3,7);
    E(3,3) = ST(2,1)^2*EQQ(1,1)+ST(1,1)^2*EQQ(4,4)+ST(2,2)^2*EQQ(2,2)...
        +ST(1,2)^2*EQQ(5,5)+ST(2,3)^2*EQQ(3,3)+ST(1,3)^2*EQQ(6,6)...
        -2*ST(2,1)*ST(1,2)*EQQ(1,5)-2*ST(1,1)*ST(2,2)*EQQ(2,4);
    E(1,2) = -ST(3,1)*ST(3,2)*EQQ(1,5)+ST(3,1)*ST(2,3)*EQQ(1,9)...
        -ST(1,1)*ST(2,1)*EQQ(7,7)...
        -ST(3,2)*ST(3,1)*EQQ(2,4)...
        -ST(1,2)*ST(2,2)*EQQ(8,8)+ST(1,2)*ST(3,3)*EQQ(6,8)...
        +ST(3,3)*ST(2,1)*EQQ(3,7)...
        +ST(1,3)*ST(3,2)*EQQ(5,9)-ST(1,3)*ST(2,3)*EQQ(9,9);
    E(1,3) = ST(2,1)*ST(3,2)*EQQ(1,5)-ST(2,1)*ST(2,3)*EQQ(1,9)...
        -ST(1,1)*ST(3,1)*EQQ(4,4)...
        +ST(2,2)*ST(3,1)*EQQ(2,4)...
        -ST(1,2)*ST(3,2)*EQQ(5,5)+ST(1,2)*ST(2,3)*EQQ(5,9)...
        -ST(2,3)*ST(2,1)*EQQ(3,7)...
        +ST(1,3)*ST(2,2)*EQQ(6,8)-ST(1,3)*ST(3,3)*EQQ(6,6);
    E(2,3) = -ST(2,1)*ST(3,1)*EQQ(1,1)+ST(2,1)*ST(1,3)*EQQ(1,9)...
        +ST(1,1)*ST(3,2)*EQQ(2,4)...
        -ST(2,2)*ST(3,2)*EQQ(2,2)...
        +ST(1,2)*ST(3,1)*EQQ(1,5)-ST(1,2)*ST(1,3)*EQQ(5,9)...
        +ST(2,3)*ST(1,1)*EQQ(3,7)-ST(2,3)*ST(3,3)*EQQ(3,3)...
        -ST(1,3)*ST(1,2)*EQQ(6,8);
    E(2,1) = E(1,2);
    E(3,1) = E(1,3);
    E(3,2) = E(2,3);
else
    E(1,1) = (ST(2,3)^2)*EQQ(5,5)+(-2*ST(2,3)*ST(3,2))*EQQ(5,9)...
        +(ST(3,2)^2)*EQQ(9,9)+(ST(1,3)^2)*EQQ(2,2)...
        +(ST(1,2)^2)*EQQ(3,3)+(ST(2,2)^2+ST(3,3)^2)*EQQ(6,6)+(-2*ST(2,2)*ST(3,3))*EQQ(6,8);
    E(2,2) = (ST(1,3)^2)*EQQ(1,1)+(-2*ST(1,3)*ST(3,1))*EQQ(1,9)...
        +(ST(3,1)^2)*EQQ(9,9)+(ST(2,3)^2)*EQQ(2,2)...
        +(ST(1,1)^2+ST(3,3)^2)*EQQ(3,3)+(-2*ST(1,1)*ST(3,3))*EQQ(3,7)+(ST(2,1)^2)*EQQ(6,6);
    E(3,3) = (ST(1,2)^2)*EQQ(1,1)+(-2*ST(1,2)*ST(2,1))*EQQ(1,5)...
        +(ST(2,1)^2)*EQQ(5,5)+(ST(1,1)^2+ST(2,2)^2)*EQQ(2,2)...
        +(-2*ST(1,1)*ST(2,2))*EQQ(2,4)+(ST(3,2)^2)*EQQ(3,3)+(ST(3,1)^2)*EQQ(6,6);
    E(1,2) = (-ST(1,3)*ST(2,3))*EQQ(1,5)+(ST(1,3)*ST(3,2))*EQQ(1,9)...
        +(ST(2,3)*ST(3,1))*EQQ(5,9)+(-ST(3,1)*ST(3,2))*EQQ(9,9)...
        +(-ST(1,3)*ST(2,3))*EQQ(2,4)+(-ST(1,1)*ST(1,2))*EQQ(3,3)...
        +(ST(1,2)*ST(3,3))*EQQ(3,7)+(-ST(2,1)*ST(2,2))*EQQ(6,6)+(ST(2,1)*ST(3,3))*EQQ(6,8);
    E(1,3) = (ST(1,2)*ST(2,3))*EQQ(1,5)+(-ST(1,2)*ST(3,2))*EQQ(1,9)...
        +(-ST(2,1)*ST(2,3))*EQQ(5,5)+(ST(2,1)*ST(3,2))*EQQ(5,9)...
        +(-ST(1,1)*ST(1,3))*EQQ(2,2)+(ST(1,3)*ST(2,2))*EQQ(2,4)...
        +(-ST(1,2)*ST(3,2))*EQQ(3,7)+(-ST(3,1)*ST(3,3))*EQQ(6,6)+(ST(2,2)*ST(3,1))*EQQ(6,8);
    E(2,3) = (-ST(1,2)*ST(1,3))*EQQ(1,1)+(ST(1,3)*ST(2,1))*EQQ(1,5)...
        +(ST(1,2)*ST(3,1))*EQQ(1,9)+(-ST(2,1)*ST(3,1))*EQQ(5,9)...
        +(-ST(2,2)*ST(2,3))*EQQ(2,2)+(ST(1,1)*ST(2,3))*EQQ(2,4)...
        +(-ST(3,2)*ST(3,3))*EQQ(3,3)+(ST(1,1)*ST(3,2))*EQQ(3,7)+(-ST(2,1)*ST(3,1))*EQQ(6,8);
    E(2,1) = E(1,2);
    E(3,1) = E(1,3);
    E(3,2) = E(2,3);
end

end


function [ E ] = EfRfRXS( EQQ, ST, s, defQS )

if defQS
    E(1,1) = -s(3)*ST(2,2)*EQQ(6,8)+s(3)*ST(3,3)*EQQ(6,6)...
        +s(2)*ST(2,2)*EQQ(6,6)-s(2)*ST(3,3)*EQQ(6,8);
    E(1,2) = s(3)*ST(2,1)*EQQ(3,7)-s(1)*ST(2,1)*EQQ(7,7);
    E(1,3) = s(2)*ST(3,1)*EQQ(2,4)-s(1)*ST(3,1)*EQQ(4,4);
    E(2,1) = s(3)*ST(1,2)*EQQ(6,8)-s(2)*ST(1,2)*EQQ(8,8);
    E(2,2) = -s(3)*ST(1,1)*EQQ(3,7)+s(3)*ST(3,3)*EQQ(3,3)...
        +s(1)*ST(1,1)*EQQ(7,7)-s(1)*ST(3,3)*EQQ(3,7);
    E(2,3) = -s(2)*ST(3,2)*EQQ(2,2)+s(1)*ST(3,2)*EQQ(2,4);
    E(3,1) = -s(3)*ST(1,3)*EQQ(6,6)+s(2)*ST(1,3)*EQQ(6,8);
    E(3,2) = -s(3)*ST(2,3)*EQQ(3,3)+s(1)*ST(2,3)*EQQ(3,7);
    E(3,3) = -s(2)*ST(1,1)*EQQ(2,4)+s(2)*ST(2,2)*EQQ(2,2)...
        +s(1)*ST(1,1)*EQQ(4,4)-s(1)*ST(2,2)*EQQ(2,4);
else
    E(1,1) = (ST(2,2)*s(2)+ST(3,3)*s(3))*EQQ(6,6) + (-ST(2,2)*s(3)-ST(3,3)*s(2))*EQQ(6,8);
    E(1,2) = (-ST(1,2)*s(1))*EQQ(3,3) + (ST(1,2)*s(3))*EQQ(3,7);
    E(1,3) = (-ST(1,3)*s(1))*EQQ(2,2) + (ST(1,3)*s(2))*EQQ(2,4);
    E(2,1) = (-ST(2,1)*s(2))*EQQ(6,6) + (ST(2,1)*s(3))*EQQ(6,8);
    E(2,2) = (ST(1,1)*s(1)+ST(3,3)*s(3))*EQQ(3,3) + (-ST(1,1)*s(3)-ST(3,3)*s(1))*EQQ(3,7);
    E(2,3) = (-ST(2,3)*s(2))*EQQ(2,2) + (ST(2,3)*s(1))*EQQ(2,4);
    E(3,1) = (-ST(3,1)*s(3))*EQQ(6,6) + (ST(3,1)*s(2))*EQQ(6,8);
    E(3,2) = (-ST(3,2)*s(3))*EQQ(3,3) + (ST(3,2)*s(1))*EQQ(3,7);
    E(3,3) = (ST(1,1)*s(1)+ST(2,2)*s(2))*EQQ(2,2) + (-ST(1,1)*s(2)-ST(2,2)*s(1))*EQQ(2,4);
end

end

