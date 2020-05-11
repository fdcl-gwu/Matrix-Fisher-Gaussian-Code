function [ EQ, EQQ, EQQQ ] = pdf_MF_moment123( s )
% calculate the first three order moments of Q

[c,dc] = pdf_MF_normal(s,1,1);

% EQ
EQ = diag(dc/c+1);

% EQQ
EQQ = zeros(3,3);
for i = 1:3
    for j = setdiff(1:3,i)
        EQQ(3*(i-1)+j,3*(i-1)+j) = EQ(i,i)*s(i)/(s(i)^2-s(j)^2)-EQ(j,j)*s(j)/(s(i)^2-s(j)^2);
        EQQ(3*(i-1)+j,3*(j-1)+i) = EQ(i,i)*s(j)/(s(i)^2-s(j)^2)-EQ(j,j)*s(i)/(s(i)^2-s(j)^2);
    end
end

EQQ(1,1) = 1-EQQ(2,2)-EQQ(3,3);
EQQ(5,5) = 1-EQQ(4,4)-EQQ(6,6);
EQQ(9,9) = 1-EQQ(7,7)-EQQ(8,8);

EQQ(1,5) = EQQ(2,4)+EQ(3,3);
EQQ(1,9) = EQQ(3,7)+EQ(2,2);
EQQ(5,9) = EQQ(6,8)+EQ(1,1);
EQQ(5,1) = EQQ(1,5);
EQQ(9,1) = EQQ(1,9);
EQQ(9,5) = EQQ(5,9);

% EQQQ
EQQQ = zeros(9,9,9);
for i = 1:3
    for j = 1:3
        for k = setdiff(1:3,j)
            if i~=j && i~=k
                EQQQ(3*(i-1)+i,3*(j-1)+k,3*(j-1)+k) = EQQ(3*(i-1)+i,3*(j-1)+j)*s(j)/(s(j)^2-s(k)^2)...
                    - EQQ(3*(i-1)+i,3*(k-1)+k)*s(k)/(s(j)^2-s(k)^2);
                EQQQ(3*(i-1)+i,3*(j-1)+k,3*(k-1)+j) = EQQ(3*(i-1)+i,3*(j-1)+j)*s(k)/(s(j)^2-s(k)^2)...
                    - EQQ(3*(i-1)+i,3*(k-1)+k)*s(j)/(s(j)^2-s(k)^2);
                
                EQQQ(3*(j-1)+k,3*(i-1)+i,3*(j-1)+k) = EQQQ(3*(i-1)+i,3*(j-1)+k,3*(j-1)+k);
                EQQQ(3*(j-1)+k,3*(i-1)+i,3*(k-1)+j) = EQQQ(3*(i-1)+i,3*(j-1)+k,3*(k-1)+j);
                
                EQQQ(3*(j-1)+k,3*(j-1)+k,3*(i-1)+i) = EQQQ(3*(i-1)+i,3*(j-1)+k,3*(j-1)+k);
                EQQQ(3*(j-1)+k,3*(k-1)+j,3*(i-1)+i) = EQQQ(3*(i-1)+i,3*(j-1)+k,3*(k-1)+j);
            end
            
            if i==j
                EQQQ(3*(j-1)+j,3*(j-1)+k,3*(j-1)+k) = EQQ(3*(j-1)+j,3*(j-1)+j)*s(j)/(s(j)^2-s(k)^2) ...
                    - EQQ(3*(j-1)+j,3*(k-1)+k)*s(k)/(s(j)^2-s(k)^2) ...
                    - EQ(j,j)*(s(j)^2+s(k)^2)/(s(j)^2-s(k)^2)^2 ...
                    + EQ(k,k)*2*s(j)*s(k)/(s(j)^2-s(k)^2)^2;
                EQQQ(3*(j-1)+j,3*(j-1)+k,3*(k-1)+j) = EQQ(3*(j-1)+j,3*(j-1)+j)*s(k)/(s(j)^2-s(k)^2) ...
                    - EQQ(3*(j-1)+j,3*(k-1)+k)*s(j)/(s(j)^2-s(k)^2) ...
                    - EQ(j,j)*2*s(j)*s(k)/(s(j)^2-s(k)^2)^2 ...
                    + EQ(k,k)*(s(j)^2+s(k)^2)/(s(j)^2-s(k)^2)^2;
                
                EQQQ(3*(j-1)+k,3*(j-1)+j,3*(j-1)+k) = EQQQ(3*(j-1)+j,3*(j-1)+k,3*(j-1)+k);
                EQQQ(3*(j-1)+k,3*(j-1)+j,3*(k-1)+j) = EQQQ(3*(j-1)+j,3*(j-1)+k,3*(k-1)+j);
                
                EQQQ(3*(j-1)+k,3*(j-1)+k,3*(j-1)+j) = EQQQ(3*(j-1)+j,3*(j-1)+k,3*(j-1)+k);
                EQQQ(3*(j-1)+k,3*(k-1)+j,3*(j-1)+j) = EQQQ(3*(j-1)+j,3*(j-1)+k,3*(k-1)+j);
            end
            
            if i==k
                EQQQ(3*(k-1)+k,3*(j-1)+k,3*(j-1)+k) = EQQ(3*(k-1)+k,3*(k-1)+k)*s(k)/(s(k)^2-s(j)^2) ...
                    - EQQ(3*(k-1)+k,3*(j-1)+j)*s(j)/(s(k)^2-s(j)^2) ...
                    - EQ(k,k)*(s(k)^2+s(j)^2)/(s(k)^2-s(j)^2)^2 ...
                    + EQ(j,j)*2*s(k)*s(j)/(s(k)^2-s(j)^2)^2;
                EQQQ(3*(k-1)+k,3*(j-1)+k,3*(k-1)+j) = EQQ(3*(k-1)+k,3*(k-1)+k)*s(j)/(s(k)^2-s(j)^2) ...
                    - EQQ(3*(k-1)+k,3*(j-1)+j)*s(k)/(s(k)^2-s(j)^2) ...
                    - EQ(k,k)*2*s(k)*s(j)/(s(k)^2-s(j)^2)^2 ...
                    + EQ(j,j)*(s(k)^2+s(j)^2)/(s(k)^2-s(j)^2)^2;
                
                EQQQ(3*(j-1)+k,3*(k-1)+k,3*(j-1)+k) = EQQQ(3*(k-1)+k,3*(j-1)+k,3*(j-1)+k);
                EQQQ(3*(j-1)+k,3*(k-1)+k,3*(k-1)+j) = EQQQ(3*(k-1)+k,3*(j-1)+k,3*(k-1)+j);
                
                EQQQ(3*(j-1)+k,3*(j-1)+k,3*(k-1)+k) = EQQQ(3*(k-1)+k,3*(j-1)+k,3*(j-1)+k);
                EQQQ(3*(j-1)+k,3*(k-1)+j,3*(k-1)+k) = EQQQ(3*(k-1)+k,3*(j-1)+k,3*(k-1)+j);
            end
        end
    end
end

for i = 1:3
    for j = setdiff([1,2,3],i)
        for k = setdiff([1,2,3],[i,j])
            EQQQ(3*(i-1)+j,3*(j-1)+k,3*(k-1)+i) = EQ(i,i)*s(j)*s(k)/(s(i)^2-s(j)^2)/(s(i)^2-s(k)^2)...
                + EQ(j,j)*s(i)*s(k)/(s(j)^2-s(i)^2)/(s(j)^2-s(k)^2)...
                + EQ(k,k)*s(i)*s(j)/(s(k)^2-s(i)^2)/(s(k)^2-s(j)^2);
            EQQQ(3*(i-1)+j,3*(j-1)+k,3*(i-1)+k) = EQ(i,i)*s(j)*s(i)/(s(i)^2-s(j)^2)/(s(i)^2-s(k)^2)...
                + EQ(j,j)*s(j)^2/(s(j)^2-s(i)^2)/(s(j)^2-s(k)^2)...
                + EQ(k,k)*s(k)*s(j)/(s(k)^2-s(i)^2)/(s(k)^2-s(j)^2);
            
            % ijkijk & ijikjk
            EQQQ(3*(i-1)+j,3*(k-1)+i,3*(j-1)+k) = EQQQ(3*(i-1)+j,3*(j-1)+k,3*(k-1)+i);
            EQQQ(3*(i-1)+j,3*(i-1)+k,3*(j-1)+k) = EQQQ(3*(i-1)+j,3*(j-1)+k,3*(i-1)+k);

            % jkijki & jkijik
            EQQQ(3*(j-1)+k,3*(i-1)+j,3*(k-1)+i) = EQQQ(3*(i-1)+j,3*(j-1)+k,3*(k-1)+i);
            EQQQ(3*(j-1)+k,3*(i-1)+j,3*(i-1)+k) = EQQQ(3*(i-1)+j,3*(j-1)+k,3*(i-1)+k);

            % jkkiij & jkikij
            EQQQ(3*(j-1)+k,3*(k-1)+i,3*(i-1)+j) = EQQQ(3*(i-1)+j,3*(j-1)+k,3*(k-1)+i);
            EQQQ(3*(j-1)+k,3*(i-1)+k,3*(i-1)+j) = EQQQ(3*(i-1)+j,3*(j-1)+k,3*(i-1)+k);

            % kiijjk & ikijjk
            EQQQ(3*(k-1)+i,3*(i-1)+j,3*(j-1)+k) = EQQQ(3*(i-1)+j,3*(j-1)+k,3*(k-1)+i);
            EQQQ(3*(i-1)+k,3*(i-1)+j,3*(j-1)+k) = EQQQ(3*(i-1)+j,3*(j-1)+k,3*(i-1)+k);

            % kijkij & ikjkij
            EQQQ(3*(k-1)+i,3*(j-1)+k,3*(i-1)+j) = EQQQ(3*(i-1)+j,3*(j-1)+k,3*(k-1)+i);
            EQQQ(3*(i-1)+k,3*(j-1)+k,3*(i-1)+j) = EQQQ(3*(i-1)+j,3*(j-1)+k,3*(i-1)+k);
        end
    end
end

EQQQ(1,1,1) = EQ(1,1)-EQQQ(1,2,2)-EQQQ(1,3,3);
EQQQ(1,1,5) = EQQQ(1,2,4)+EQQ(1,9);
EQQQ(1,1,9) = EQQQ(1,3,7)+EQQ(1,5);
EQQQ(1,5,5) = EQ(1,1)-EQQQ(1,4,4)-EQQQ(1,6,6);
EQQQ(1,5,9) = EQQQ(1,6,8)+EQQ(1,1);
EQQQ(1,9,9) = EQ(1,1)-EQQQ(1,7,7)-EQQQ(1,8,8);
EQQQ(5,5,5) = EQ(2,2)-EQQQ(5,4,4)-EQQQ(5,6,6);
EQQQ(5,5,9) = EQQQ(5,6,8)+EQQ(5,1);
EQQQ(5,9,9) = EQ(2,2)-EQQQ(5,7,7)-EQQQ(5,8,8);
EQQQ(9,9,9) = EQ(3,3)-EQQQ(9,7,7)-EQQQ(9,8,8);

EQQQ(1,5,1) = EQQQ(1,1,5);
EQQQ(5,1,1) = EQQQ(1,1,5);

EQQQ(1,9,1) = EQQQ(1,1,9);
EQQQ(9,1,1) = EQQQ(1,1,9);

EQQQ(5,1,5) = EQQQ(1,5,5);
EQQQ(5,5,1) = EQQQ(1,5,5);

EQQQ(1,9,5) = EQQQ(1,5,9);
EQQQ(5,1,9) = EQQQ(1,5,9);
EQQQ(5,9,1) = EQQQ(1,5,9);
EQQQ(9,1,5) = EQQQ(1,5,9);
EQQQ(9,5,1) = EQQQ(1,5,9);

EQQQ(9,1,9) = EQQQ(1,9,9);
EQQQ(9,9,1) = EQQQ(1,9,9);

EQQQ(5,9,5) = EQQQ(5,5,9);
EQQQ(9,5,5) = EQQQ(5,5,9);

EQQQ(9,5,9) = EQQQ(5,9,9);
EQQQ(9,9,5) = EQQQ(5,9,9);

end

