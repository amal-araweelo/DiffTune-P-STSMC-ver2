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

function ud = controller(X, Xref, k_vec, theta_r_dot, theta_r_2dot, J_m, N, dt) % t for time
    global v;

    % Controller gains
    k1 = k_vec(1);
    k2 = k_vec(2);
    k_pos = k_vec(3);

    % States
    omega_m = X(1);
    omega_l = X(2);
    theta_l = X(4);
    theta_r = Xref(4);
    
    % Controllers
 
    % P-controller
    omega_r = k_pos * (theta_r - theta_l) + N * theta_r_dot;
    omega_r_dot = k_pos * (theta_r_dot - omega_l) + N * theta_r_2dot;
   
    % PI-controller
    % u = int((omega_r - omega_m) * k_vel * 1/tau_i) + (omega_r - omega_m) * k_vel + diff(omega_r) * J_m;

    % STSMC controller
    s = omega_m - omega_r; % Error
    v_dot = -k2 * sgn_approx(100*s);
    
    if (isempty(v)) % initialise v to zero in first iteration
        v = 0;
    end
    v = v + v_dot * dt;

    u_smc = -k1 * sqrt(abs(s)) * sgn_approx(s) + v;
    % fprintf('omega_r_dot = %d \n', omega_r_dot);
    u = u_smc + J_m * omega_r_dot;
    % fprintf('u = %d \n', u);

    % Output
    ud = u;
    % fprintf('ud = %d \n', ud);
end 