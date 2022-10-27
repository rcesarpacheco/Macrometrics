function C = recover_level_C(delta_c,delta_hat,C0)
C = delta_c;
for i=1:length(delta_c)
    if i==1
        C(i) = C0+delta_c(i)+delta_hat;
    else
        C(i) = C(i-1)+delta_c(i)+delta_hat;
    end
end