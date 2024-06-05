% define the controller here

% Inputs:
%   theta_r: load position reference
%   omega_r: motor velocity reference
%   omega_dot_r
%   u: torque command

% States and Outputs
%   omega_m: motor angular velocity
%   theta_m: motor angular position
%   omega_l: load angular velocity
%   theta_l: load angular position
%   ud: containing motor/load angular velocity/position (4outputs)

% X(1) = omega_m
% X(2) = theta_m
% X(3) = omega_l
% X(4) = theta_l

% Xref(1) = theta_r
% Xref(2) = omega_r

function ud = controller(X, Xref, k_vec)

% Controller gains
k1 = k_vec(1);
k2 = k_vec(2);
k_pos = k_vec(3);


% Controllers

% P-controller
omega_r = k_pos * (Xref(1) - X(4)) + N * diff(theta_r);

% PI-controller
% u = int((omega_r - omega_m) * k_vel * 1/tau_i) + (omega_r - omega_m) * k_vel + diff(omega_r) * J_m;

% STSMC controller
s = X(1) - omega_r;
sgn(x) = 2/pi*atan(10*x);
v_dot = -k2 * sgn(s);
u_smc = -k1 * sqrt(abs(s)) * sgn(s) + int(v_dot);
u = u_smc + J_m * diff(omega_r);

ud = u;
end