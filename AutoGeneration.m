% This script uses CasADi to autogenerate functions for online Jacobian
% evalutaion

clear all;
clc;
import casadi.*;

%% Define the dimensions
dim_state = 4; % dimension of system state (omega_m, theta_m, omega_r, theta_r)
dim_control = 1;  % dimension of control inputs (u, theta_r, omega_r)
dim_controllerParameters = 3;  % dimension of controller parameters (k_1, k_2, k_pos)

%% Load constant physical parameters
% Sampling time
dt = MX.sym('dt', 1); % (should be set to 1-8 kHz in runDiffTune.m)
t = MX.sym('t', 1);
v_temp = MX.sym('v_temp', 1);

% Constant drive train parameters
N = MX.sym('N', 1);             % N: Gearing ratio
J_m = MX.sym('J_m', 1);         % J_m: Motor inertia
J_l = MX.sym('J_l', 1);         % J_l: Load inertia
K_S = MX.sym('K_S', 1);         % K_s: Shaft stifness
D_S = MX.sym('D_S', 1);         % D_s: Shaft damping coefficinet
T_Cm = MX.sym('T_Cm', 1);       % T_Cm: Motor Coulomb friction
T_Cl = MX.sym('T_Cl', 1);
beta_m = MX.sym('beta_m', 1);   % beta_m: Motor viscous friction coefficient
beta_l = MX.sym('beta_l', 1);

param = [N J_m J_l K_S D_S T_Cm T_Cl beta_m beta_l];

%% casADI-lize all the variables in the computation
X = MX.sym('X', dim_state);        % system state
Xref = MX.sym('Xref', 1);  % system reference state

% Desired values
theta_r_dot = MX.sym('theta_r_dot', 1);
theta_r_2dot = MX.sym('theta_r_2dot', 1);

%% k is the collection of controller parameters 
k_vec = MX.sym('k_vec', dim_controllerParameters); % gains for P-STSMC
%% Define the control input
u = MX.sym('u', dim_control);

%% Define the dynamics (discretized via Forward Euler)
dynamics = X + dt * dynamics(t, X, u, param);
                    
%% Compute the control action, denoted by h
[h, v] = controller(X, Xref, k_vec, theta_r_dot, theta_r_2dot, param, dt, v_temp); % Xref(4) = theta_r is the desired trajectory

%% Generate jacobians
grad_f_X = jacobian(dynamics,X);
grad_f_u = jacobian(dynamics,u);
grad_h_X = jacobian(h,X);
grad_h_theta = jacobian(h,k_vec);

%% Function-lize the generated jacobians

% inputs_f denotes the input arguments to the dynamics
grad_f_X_fcn = Function('grad_f_X_fcn',{X, dt, u, N, J_m, J_l, K_S, D_S, T_Cm, T_Cl, beta_m, beta_l},{grad_f_X});
grad_f_u_fcn = Function('grad_f_u_fcn',{X, dt, u, N, J_m, J_l, K_S, D_S, T_Cm, T_Cl, beta_m, beta_l},{grad_f_u});

% inputs_h denote the input arguments to the controller
grad_h_X_fcn = Function('grad_h_X_fcn',{X, Xref, k_vec, theta_r_dot, theta_r_2dot, N, J_m, dt, v_temp},{grad_h_X});
grad_h_theta_fcn = Function('grad_h_theta_fcn',{X, Xref, k_vec, theta_r_dot, theta_r_2dot, N, J_m, dt, v_temp},{grad_h_theta});

%% Generate mex functions
opts = struct('main', true,...
              'mex', true);

mkdir mex
cd('./mex');
grad_f_X_fcn.generate('grad_f_X_fcn.c',opts);
grad_f_u_fcn.generate('grad_f_u_fcn.c',opts);
grad_h_X_fcn.generate('grad_h_X_fcn.c',opts);
grad_h_theta_fcn.generate('grad_h_theta_fcn.c',opts);

mex grad_f_X_fcn.c -largeArrayDims
mex grad_f_u_fcn.c -largeArrayDims
mex grad_h_X_fcn.c -largeArrayDims
mex grad_h_theta_fcn.c -largeArrayDims
cd('..\')