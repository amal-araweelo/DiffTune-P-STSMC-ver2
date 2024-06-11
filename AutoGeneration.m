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
dt = MX.sym('dt',1); % (should be set to 1-8 kHz in runDiffTune.m)
t = MX.sym('t', 1);

% Constant drive train parameters
N = MX.sym('N',1);              % N: Gearing ratio
J_m = MX.sym('J_m', 1);         % J_m: Motor inertia
J_l = MX.sym('J_l', 1);         % J_l: Load inertia
K_S = MX.sym('K_S', 1);         % K_s: Shaft stifness
D_S = MX.sym('D_S', 1);         % D_s: Shaft damping coefficinet
T_C = MX.sym('T_C', 1);       % T_Cm: Motor Coulomb friction
% T_S = MX.sym('T_S', 1);       % T_Sm: Motor static friction coefficient
b_fr = MX.sym('b_fr', 1);   % beta_m: Motor viscous friction coefficient

% Disturbances
% T_Fm = MX.sym('T_Fm', 1);       % T_Fm: Motor friction
% T_Fl = MX.sym('T_Fl', 1);       % T_Fl: Load friction
% T_l = MX.sym('T_l', 1);         % T_l: Load torque

param = [N J_m J_l K_S D_S T_C b_fr];

%% casADI-lize all the variables in the computation
X = MX.sym('X', dim_state);        % system state
Xref = MX.sym('Xref', 1);  % system reference state

% Elementwise split of the necessary states
% omega_m = X(1);
% omega_l = X(2);

% Desired values
theta_r_dot = MX.sym('theta_r_dot', 1);
theta_r_2dot = MX.sym('theta_r_2dot', 1);

%% k is the collection of controller parameters 
k_vec = MX.sym('k_vec', dim_controllerParameters); % gains for P-STSMC

% Split into elementwise control parameters
% k1 = k_vec(1);
% k2 = k_vec(2);
% k_pos = k_vec(3);


%% Define the control input
u = MX.sym('u', dim_control);

%% Define the dynamics (discretized via Forward Euler)
dynamics = X + dt * dynamics(t, X, u, param);
                    
%% Compute the control action, denoted by h
h = controller(X, Xref, k_vec, theta_r_dot, theta_r_2dot, param, dt); % Xref(4) = theta_r is the desired trajectory

%% Generate jacobians
grad_f_X = jacobian(dynamics,X);
grad_f_u = jacobian(dynamics,u);
grad_h_X = jacobian(h,X);
grad_h_theta = jacobian(h,k_vec);

%% Function-lize the generated jacobians

% inputs_f denotes the input arguments to the dynamics
grad_f_X_fcn = Function('grad_f_X_fcn',{X, dt, u, N, J_m, J_l, K_S, D_S, T_C, b_fr},{grad_f_X});
grad_f_u_fcn = Function('grad_f_u_fcn',{X, dt, u, N, J_m, J_l, K_S, D_S, T_C, b_fr},{grad_f_u});

% inputs_h denote the input arguments to the controller
grad_h_X_fcn = Function('grad_h_X_fcn',{X, Xref, k_vec, theta_r_dot, theta_r_2dot, J_m, N, dt},{grad_h_X});
grad_h_theta_fcn = Function('grad_h_theta_fcn',{X, Xref, k_vec, theta_r_dot, theta_r_2dot, J_m, N, dt},{grad_h_theta});

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