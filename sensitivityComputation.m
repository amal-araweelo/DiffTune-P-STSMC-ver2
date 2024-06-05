%% Sensivity computations

% This function defines how the sensitivity propgation runs at each sample time

% sensitivity: dxdtheta
% X: actual state (states: omega_m, omega_l, theta_m, theta_l)
% Xref: desired/reference state
% u: torque command (output of controller, input to system)
% J_m: moment of inertia (from motor mechanical parameters)
% k: controller gains (k_vec = [k1;k2;k_pos])
% dt: taken from simulation parameter
% dXdphi and dudphi: calculated sensitivities

function [dXdphi,dudphi] = sensitivityComputation(sensitivity,X,Xref,u,J_m,k,dt)

% evaluate the Jacobians
dfdX = grad_f_X_fcn(X,u,dt,J_m);
dfdX = full(dfdX);    % full() converts sparse matrix to full matrix (see https://se.mathworks.com/help/matlab/ref/full.html)

dfdu = grad_f_u_fcn(X,u,dt,J_m);
dfdu = full(dfdu);

dhdX = grad_h_X_fcn(X,Xref,k,J_m);
dhdX = full(dhdX);

dhdtheta = grad_h_theta_fcn(X,Xref,k,J_m);
dhdtheta = full(dhdtheta);

% assemble the Jacobians to compute the sensitivity
dXdphi = (dfdX + dfdu * dhdX) * sensitivity + dfdu * dhdtheta;
dudphi = dhdX * sensitivity + dhdtheta;

end





