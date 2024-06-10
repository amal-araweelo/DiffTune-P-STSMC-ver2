% Define the system dynamics (continuous time) here. 
% It will be used in propagating the continuous dynamics.
% X = [omega_m; omega_l; theta_m; theta_l]

function dXdt = dynamics(t,X,u,param)
dXdt = zeros(4,1);  % where 4 comes from dim_state

omega_m = X(1);
omega_l = X(2);
theta_m = X(3);
theta_l = X(4);

N = param.N;
J_m = param.J_m;
J_l = param.J_l;
K_S = param.K_S;
D_S = param.D_S;
b_fr = param.b_fr;
T_C = param.T_C;

T_l = K_S*(theta_m/N - theta_l) + D_S*(omega_m/N - omega_l);
T_Fm = omega_m*b_fr + sgn_approx(omega_m*10)*T_C;
T_Fl = omega_l*b_fr + sgn_approx(omega_l*10)*T_C + 0;

dXdt(1) = 1/J_m*u - 1/J_m*T_Fm - 1/(N*J_m)*T_l; % omega_m_dot
dXdt(2) = omega_m;                              % theta_m_dot
dXdt(3) = T_l/J_l - T_Fl/J_l;                   % omega_l_dot
dXdt(4) = omega_l;                              % theta_l_dot

end