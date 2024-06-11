% Define the system dynamics (continuous time) here. 
% It will be used in propagating the continuous dynamics.
% X = [omega_m; omega_l; theta_m; theta_l]

function dXdt = dynamics(t, X, u, param)
dXdt = zeros(4, 1); % where (4, 1) comes from (dim_state, dim_control)

% States (X = [omega_m; omega_l; theta_m; theta_l])
omega_m = X(1);
omega_l = X(2);
theta_m = X(3);
theta_l = X(4);

% Parameters (param = [N J_m J_l K_S D_S T_C T_S b_fr])
N = param(1);
J_m = param(2);
J_l = param(3);
K_S = param(4);
D_S = param(5);
T_C = param(6);
% T_S = param(7);
b_fr = param(8); 

% Disturbances computations
T_l = K_S * (theta_m / N - theta_l) + D_S * (omega_m / N - omega_l);
T_Fm = omega_m * b_fr + sgn_approx(omega_m * 10) * T_C;
T_Fl = omega_l * b_fr + sgn_approx(omega_l * 10) * T_C + 0;

% Derivative computations
dXdt(1) = 1/J_m*u - 1/J_m*T_Fm - 1/(N*J_m)*T_l; % omega_m_dot
dXdt(2) = T_l/J_l - T_Fl/J_l;                   % omega_l_dot
dXdt(3) = omega_m;                              % theta_m_dot
dXdt(4) = omega_l;                              % theta_l_dot

end