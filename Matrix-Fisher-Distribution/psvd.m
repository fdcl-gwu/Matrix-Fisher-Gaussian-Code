function [ U, S, V ]=psvd( F )
%PSVD: Proper singular value decomposition
%   [U S V]=psvd(F) is the proper singular value decomposition of a 3x3
%   matrix F. It returns two rotation matrices U, V in SO(3), and a
%   diagonal matrix S of the proper singular values such that F=U*S*V'
%
%   See T. Lee, "Bayesian Attitude Estimation with the Matrix Fisher
%   Distribution on SO(3)", 2017, http://arxiv.org/abs/1710.03746

[U,S,V]=svd(F);
for i = 1:3
    temp = nonzeros(U(:,i));
    if temp(1) < 0
        U(:,i) = -U(:,i);
        V(:,i) = -V(:,i);
    end
end

S=S*diag([1 1 det(U*V)]);
U=U*diag([1 1 det(U)]);
V=V*diag([1 1 det(V)]);

end

