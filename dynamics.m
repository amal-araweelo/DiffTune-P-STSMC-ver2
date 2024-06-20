% Define the system dynamics (continuous time) here. 
% It will be used in propagating the continuous dynamics.
% X = [omega_m; omega_l; theta_m; theta_l]

function dXdt = dynamics(t, X, u, param)

% States
omega_m = X(1);
omega_l = X(2);
theta_m = X(3);
theta_l = X(4);

% Parameters (param = [N J_m J_l K_S D_S T_Cm T_Cl beta_m beta_l])
N = param(1);
J_m = param(2);
J_l = param(3);
K_S = param(4);
D_S = param(5);
T_Cm = param(6);
T_Cl = param(7);
beta_m = param(8);
beta_l = param(9);

% Disturbances computations
T_l = K_S * (1/N * theta_m - theta_l) + D_S * (1/N * omega_m - omega_l);
T_Fm = omega_m * beta_m + sgn_approx(100 * omega_m) * T_Cm;
T_Fl = omega_l * beta_l + sgn_approx(100 * omega_l) * T_Cl;

% Derivative computations
% dXdt1 = 1/J_m*u - 1/J_m*T_Fm - 1/(N*J_m)*T_l; % omega_m_dot
% dXdt2 = T_l/J_l - T_Fl/J_l;                   % omega_l_dot
% dXdt3 = omega_m;                              % theta_m_dot
% dXdt4 = omega_l;                              % theta_l_dot

dXdt1 = 1/J_m * (u - T_Fm - 1/N * T_l); % omega_m_dot
dXdt2 = 1/J_l * (T_l - T_Fl);           % omega_l_dot
dXdt3 = omega_m;                        % theta_m_dot
dXdt4 = omega_l;                        % theta_l_dot

dXdt = [dXdt1; dXdt2; dXdt3; dXdt4];

end