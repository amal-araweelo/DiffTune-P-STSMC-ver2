% define the controller here

% Inputs:
%   theta_r: load position reference
%   omega_r: motor velocity reference
%   omega_dot_r
%   u: torque command

% States include
% omega_m: Motor angular velocity
% omega_l: Load angular velocity
% theta_m: motor angular position
% theta_l: load angular position
% X = [omega_m; omega_l; theta_m; theta_l]
% Xref (ini) = [omega_m; omega_l; theta_m; theta_r]

% ud: containing motor/load angular velocity/position (4outputs)

function [ud, v] = controller(X, Xref, k_vec, theta_r_dot, theta_r_2dot, param, dt, v) % t for time

% Controller gains
k1 = k_vec(1);
k2 = k_vec(2);
k_pos = k_vec(3);

% States
omega_m = X(1);
omega_l = X(2);
theta_l = X(4);
theta_r = Xref;

% Parameters (param = [N J_m J_l K_S D_S T_Cm T_Cl beta_m beta_l])
N = param(1);
J_m = param(2);

% P-controller
omega_r = k_pos * (theta_r - theta_l) + N * theta_r_dot;
omega_r_dot = k_pos * (theta_r_dot - omega_l) + N * theta_r_2dot;
% omega_r = k_pos * (theta_l - theta_r) + theta_r_dot;            % (Table 5.2)
% omega_r_dot = k_pos * (omega_l - theta_r_dot) + theta_r_2dot;   % (5.9)

% STSMC controller
s = omega_m - omega_r; % Error
v_dot = - k2 * sgn_approx(100*s);

v = v + v_dot * dt;

u_smc = - k1 * sqrt(abs(s)) * sgn_approx(100*s) + v;
u = u_smc + J_m * omega_r_dot;

% Output
ud = u;

end 