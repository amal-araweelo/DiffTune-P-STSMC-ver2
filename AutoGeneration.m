% This script uses CasADi to autogenerate functions for online Jacobian
% evalutaion

clear all;
import casadi.*;

%% Define the dimensions
dim_state = 4; % dimension of system state (omega_m, theta_m, omega_r, theta_r)
dim_control = 3;  % dimension of control inputs (u, theta_r, omega_r)
dim_controllerParameters = 3;  % dimension of controller parameters (k_1, k_2, k_pos)

%% Load constant physical parameters
% Sampling time
dt = MX.sym('dt',1); % (should be set to 1-8 kHz in runDiffTune.m)

% Constant drive train parameters
N = MX.sym('M',1);              % N: Gearing ratio
J_m = MX.sym('J_m', 1);         % J_m: Motor inertia
J_l = MX.sym('J_l', 1);         % J_l: Load inertia
% K_s = MX.sym('K_s', 1);         % K_s: Shaft stifness
% D_s = MX.sym('D_s', 1);         % D_s: Shaft damping coefficinet
% T_Cm = MX.sym('T_Cm', 1);       % T_Cm: Motor Coulomb friction
% T_Sm = MX.sym('T_Sm', 1);       % T_Sm: Motor static friction coefficient
% omega_s = MX.sym('omega_s', 1); % omega_s: Motor Stribeck velocity
% beta_m = MX.sym('beta_m', 1);   % beta_m: Motor viscous friction coefficient

% Disturbances
%d_e = MX.sym('d_e', 1);         % d_e: Input torque ripples and harmonics
T_Fm = MX.sym('T_Fm', 1);       % T_Fm: Motor friction
T_Fl = MX.sym('T_Fl', 1);       % T_Fl: Load friction
T_l = MX.sym('T_l', 1);         % T_l: Load torque

%% casADI-lize all the variables in the computation
X = MX.sym('X',dim_state);          % system state
Xref = MX.sym('X_ref', dim_state);  % system reference state

% Elementwise split of the necessary states
omega_m = X(1);
omega_l = X(2);
theta_r = Xref(4);

% Load the desired values into a struct
desired = MX.sym('theta_r', 1);
theta_r_dot = MX.sym('theta_r_dot', 1);

%% k is the collection of controller parameters 
k_vec = MX.sym('k_vec',dim_controllerParameters); % gains for P-STSMC

% Split into elementwise control parameters
k1 = k_vec(1);
k2 = k_vec(2);
k_pos = k_vec(3);


%% Define the control input
% u = MX.sym('u',dim_control);
u = MX.sym('u', 1);    % ændret fordi det u vi bruger her er inputtet til systemet u og er 1 dimensionelt.

%% Define the dynamics (discretized via Forward Euler)

dynamics = X + dt * [1/J_m*u - 1/J_m*T_Fm - 1/(N*J_m)*T_l;
                    omega_m;
                    T_l/J_l - T_Fl/J_l;
                    omega_l]; 
                    
%% Compute the control action, denoted by h
h = controller(X, Xref, k_vec, theta_r_dot, J_m, N, dt); % Xref(4) = theta_r is the desired trajectory

%% Generate jacobians
grad_f_X = jacobian(dynamics,X);
grad_f_u = jacobian(dynamics,u);
grad_h_X = jacobian(h,X);
grad_h_theta = jacobian(h,k_vec);

%% Function-lize the generated jacobians

% inputs_f denotes the input arguments to the dynamics
grad_f_X_fcn = Function('grad_f_X_fcn',{X, dt, u, J_m, N, J_l},{grad_f_X});
grad_f_u_fcn = Function('grad_h_u_fcn',{X, dt, u, J_m, N, J_l},{grad_f_u});

% inputs_h denote the input arguments to the controller
grad_h_X_fcn = Function('grad_h_X_fcn',{X, Xref, k_vec, theta_r_dot, J_m, dt},{grad_h_X});
grad_h_theta_fcn = Function('grad_h_theta_fcn',{X, Xref, k_vec, theta_r_dot, J_m, dt},{grad_h_theta});

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