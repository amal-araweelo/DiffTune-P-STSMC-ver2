% Define the system dynamics (continuous time) here. 
% It will be used in propagating the continuous dynamics.
% X = [omega_m; omega_l; theta_m; theta_l]

function dXdt = dynamics(t,X,u,param)
omega_m = X(1);
omega_l = X(2);
theta_m = X(3);
theta_l = X(4);

N = param.N;
J_m = param.J_m;
J_l = param.J_l;

T_l = param.K_S*(theta_m/param.N - theta_l) + param.D_S*(omega_m/N - omega_l);
T_Fm = omega_m*param.b_fr + sgn_approx(omega_m*10)*param.T_C;
T_Fl = omega_l*param.b_fr + sgn_approx(omega_l*10)*param.T_C + 0;

dXdt(1) = 1/J_m*u - 1/J_m*T_Fm - 1/(N*J_m)*T_l; % omega_m_dot
dXdt(2) = omega_m;                              % theta_m_dot
dXdt(3) = T_l/J_l - T_Fl/J_l;                   % omega_l_dot
dXdt(4) = omega_l;                              % theta_l_dot

end