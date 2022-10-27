function [Z,T,Q] = state_space_representation(par)
% Kalman Filter: conventions
% a_t=T*a_{t-1}+eta_t. eta_t~N(0,Q)
% yt=Z*a_t+e_t e_t~N(0,H)
gamma1 = par(1);
gamma2 = par(2);
gamma3 = par(3);
gamma4 = par(4);

psi11 = par(5);
psi12 = par(6);
psi13 = par(7);
psi14 = par(8);

psi21 = par(9);
psi22 = par(10);
psi23 = par(11);
psi24 = par(12);

phi1 = par(13);
phi2 = par(14);

var1 = par(15);
var2 = par(16);
var3 = par(17);
var4 = par(18);
Z = [gamma1,0,1,0,0,0,0,0,0,0;
    gamma2,0,0,0,1,0,0,0,0,0;
    gamma3,0,0,0,0,0,1,0,0,0;
    gamma4,0,0,0,0,0,0,0,1,0];

T = [phi1,phi2,0,0,0,0,0,0,0,0;
    1,0,0,0,0,0,0,0,0,0;
    0,0,psi11,psi21,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,0,psi12,psi22,0,0,0,0;
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,0,psi13,psi23,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,0,psi14,psi24;
    0,0,0,0,0,0,0,0,1,0];
Q = diag([1,0,var1,0,var2,0,var3,0,var4,0]);
end
