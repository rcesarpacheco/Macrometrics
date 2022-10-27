function ll = log_likelihood(a0,P0,Z,T,Q,H,y)
kf = kalman_filter(a0,P0,Z,T,Q,H,y);
inov = kf.inov;
cov_inov_matrix = kf.cov_inov_matrix;
nt = size(y,2);
ll=-log(2*pi)*nt*size(y,1)/2;
for i=1:nt
    ll = ll - 1/2*log(det(cov_inov_matrix(:,:,i)))-1/2*inov(:,i)'*inv(cov_inov_matrix(:,:,i))*inov(:,i);
end



