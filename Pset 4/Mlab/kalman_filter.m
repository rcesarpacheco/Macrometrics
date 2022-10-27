% function [at,K,inov,cov_inov_matrix] = kalman_filter(a0,P0,Z,T,Q,H,y)
function kf = kalman_filter(a0,P0,Z,T,Q,H,y)
nt = size(y,2);
at = zeros(size(T,1),nt);
P = zeros(size(Q,1),size(Q,2),nt);
K = zeros(size(P,1),size(y,1),nt);
cov_inov_matrix = zeros(size(y,1),size(y,1),nt);
for i=1:nt
    if i==1
        atm = T*a0;
        ptm = T*P0*T'+Q;
    else
        atm = T*at(:,i-1);
        ptm = T*P(:,:,i-1)*T'+Q;
    end
        F = Z*ptm*Z'+H;
        ytm = Z*atm;
        inov(:,i) = y(:,i)-ytm;
        K(:,:,i) = ptm*Z'*inv(F);
        P(:,:,i)=ptm-K(:,:,i)*Z*ptm;
        at(:,i)=atm+K(:,:,i)*inov(:,i);
        cov_inov_matrix(:,:,i) = F;
end
kf.at = at;
kf.K = K;
kf.inov = inov;
kf.cov_inov_matrix = cov_inov_matrix;
end
