function obj = obj_function_ml(par,a0,y)
[Z,T,Q] = state_space_representation(par);
P0 = T*Q*T';
H= zeros(size(y,1),size(y,1));
obj = -log_likelihood(a0,P0,Z,T,Q,H,y);
end