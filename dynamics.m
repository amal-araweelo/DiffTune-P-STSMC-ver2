% Define the system dynamics (continuous time) here. 
% It will be used in propagating the continuous dynamics.

function dXdt = dynamics(t,X,u,param)

J_m = param.J_m;
T_Fm = param.T_Fm;
T_l = param.T_l;
J_l = param.J_L;
T_Fl = param.F_fl;

dXdt(1) = 1/J_m*u - 1/J_m*T_Fm - 1/(N*J_m)*T_l; % omega_m_dot
dXdt(2) = omega_m;                              % theta_m_dot
dXdt(3) = T_l/J_l - T_Fl/J_l;                   % omega_l_dot
dXdt(4) = omega_l;                              % theta_l_dot

end