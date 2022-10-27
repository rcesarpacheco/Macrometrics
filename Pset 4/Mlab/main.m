% Originals paramaters of the model
gamma1 = 0.717;
gamma2 = 0.521;
gamma3 = 0.47;
gamma4 = 0.602;

psi11 = -0.04;
psi12 = -0.087;
psi13 = -0.414;
psi14 = 0.108;

psi21 = -0.137;
psi22 = 0.154;
psi23 = -0.206;
psi24 = 0.448;

phi1 = 0.545;
phi2 = 0.032;
var1 = (0.488*1e-2)^2;
var2 = (0.769*1e-2)^2;
var3 = (0.735*1e-2)^2;
var4 = (0.540*1e-2)^2;


par = [gamma1, gamma2, gamma3, gamma4,...
    psi11,psi12,psi13,psi14,...
    psi21,psi22,psi23,psi24,...
    phi1,phi2,var1,var2,var3,var4];
%% part b
n_states=10;
init_par = par;
a0 = zeros(n_states,1);
options = optimset('MaxFunEvals',1e5,'MaxIter',1e5,'Display','iter');
new_y = readmatrix("y_fred.csv");
params_estimated = fminsearch(@(x) obj_function_ml(x,a0,new_y),init_par,options)
%% export results to table
params_estimated(15:18)=sqrt(par(15:18));
var_names = ["gamma1","gamma2","gamma3","gamma4","psi11","psi12","psi13","psi14","psi21","psi22","psi23","psi24","phi1","phi2","sigma1","sigma2","sigma3","sigma4"]
table_results = table(var_names',params_estimated')

results_matrix = [params_estimated(1:4);params_estimated(5:8); params_estimated(9:12);params_estimated(15:18)]
csvwrite('params_estimated_2b.csv',results_matrix,m)