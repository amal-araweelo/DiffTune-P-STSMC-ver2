% Use this script to 
% run the simulation with DiffTune

% Constant drive train parameters
% N: Gearing ratio
% J_m: Motor inertia
% J_l: Load inertia
% K_s: Shaft stifness
% D_s: Shaft damping coefficinet
% T_Cm: Motor Coulomb friction
% T_Sm: Motor static friction coefficient
% omega_s: Motor Stribeck velocity
% beta_m: Motor viscous friction coefficient

% Disturbances
% d_e: Input torque ripples and harmonics
% T_Fm: Motor friction
% T_Fl: Load friction
% T_l: Load torque

% States include
% omega_m: Motor angular velocity
% omega_l: Load angular velocity
% theta_m: motor angular position
% theta_l: load angular position

% Control includes
% Theta_r: Load postition reference
% Omega_r: Motor velocity reference
% u: Torque command

close all;
clear all;

addpath('mex\');
import casadi.*

%% define the dimensions
dim_state = 4; % dimension of system state
dim_control = 1;  % dimension of control inputs
dim_controllerParameters = 3??;  % dimension of controller parameters

%% Define simulation parameters (e.g., sample time dt, duration, etc)
%dt = 0.01; %0.125ms
%time = 0:dt:10;

%% constant parameters
% omega_s: Motor Stribeck velocity
% Motor electrical parameters
r_s = 3.6644;           % Ohm -- Stator winding resistance (per phase)
L_d = 21.4e-3;          % mH -- Rotating field inductance
L_q = 1.2*L_d;          % mH -- Rotating torque inductance
P = 6;                  % Non-dimensional -- Number of poles
k_T = 1.43;             % Nm/A -- Torque constant
k_E = 87*2*pi*60*0.001; % Vs/rad -- Voltage constant
lambda_m = 0.3148;      % Vs/rad -- amplitude of the flux linkages
                        %   established by the permanent magnet as viewed
                        %   from the stator phase windings.
% Motor mechanical parameters
J_m = 2.81e-4 + 5.5e-4; % kgm^2 -- Moment of inertia
N = 1;                  % -- Gear ratio
% Values of friction and shaft parameters
% Taken from Table 4.3: Summary of calculated friction and shaft parameters
% (page 40, Dimitrios Papageorgiou phd thesis)
% Shaft constants
K_S = 32.94;    % N m rad^(-1)
D_S = 0.0548;   % N m s rad^(-1)
% Coulomb friction
% (assuming T_C is the average of T_C_m and T_C_l)
T_C = (0.0223 + 0.0232) / 2;    % N m
% Static friction
% (assuming T_S is the average of T_S_m and T_S_l)
T_S = (0.0441 + 0.0453) / 2;    % N m
% Friction constants
b_fr = 0.0016;  % N m s rad^(-1)
% Load inertia      (not sure...)
J_l = 1; % kgm^2 -- Moment of inertia
inv_J_l = J_l;


%% Initialize controller gains (must be a vector of size dim_controllerParameters x 1)
% STSMC (in nonlinear controller for omega_m)
k1 = 1.453488372 * 2.45 * 0.99; % use proportional gain from PI controller (k_vel = 1.45*2.45)
k2 = 50;
k_pos = 25;      % ignored when hand-tuning STSMC
k_vec = [k1, k2, k_pos];


%% Define desired trajectory if necessary
theta_r = sin(2*pi*time);   % theta_r is a sine wave with frequency 1 Hz
theta_r_dot = 2 * pi * cos(2*pi*time);


% skriver af fra ligningerne der bestemmer dem
T_l = K_s*(theta_m/N - theta_l) + D_s*(omega_m/N - omega_l);
T_Fm = omega_m*b_fr + sgn(omega_m*10)*T_C;
temp = abs(atan(omega_l))*pi/2;
    if temp > 1
        temp = 1;
    elseif temp < 0
        temp = 0;
    end
T_Fl = omega_l*b_fr + sgn(omega_l*10)*T_C + 0;


%% Initialize variables for DiffTune iterations
learningRate = 2;  
maxIterations = 100;
itr = 0;

loss_hist = [];  % storage of the loss value in each iteration
% rmse_hist = []; % If we want video
param_hist = []; % storage of the parameter value in each iteration
gradientUpdate = zeros(dim_controllerParameters,1); % define the parameter update at each iteration

%% DiffTune iterations
%while (1)
    itr = itr + 1;

    % initialize state
    X_storage = zeros(dim_state,1);
    
    % initialize sensitivity
    dx_dtheta = zeros(dim_state,dim_controllerParameters);
    du_dtheta = zeros(dim_control,dim_controllerParameters);

    % initializegradient of loss
    theta_gradient = zeros(1,dim_controllerParameters);

    % initialize reference state if necessary

    for k = 1:length(time)-1
        % load current state and current reference
        X = X_storage(:,end);
        Xref = Xref_storage(:,end);
 
        % compute the control action (of dimension dim_control x 1)
        u = controller(...);

        compute the sensitivity 
        %[dx_dtheta, du_dtheta] = sensitivityComputation(dx_dtheta,...);
        
        % accumulating the gradient of loss wrt controller parameters
        % you need to provide dloss_dx and dloss_du here
        theta_gradient = theta_gradient + dloss_dx * dx_dtheta + + dloss_du * du_dtheta;

        % integrate the ode dynamics
        % [~,sold] = ode45(@(t,X) dynamics(t,...),[time(k) time(k+1)],X);
        % X_storage = [X_storage sold(end,:)'];

        % integrate the reference system if necessary
        
    % end
    
    % loss is the squared norm of the position tracking error
    % loss = ...
    % loss_hist = [loss_hist loss];

    % update the gradient
    % gradientUpdate = - learningRate * theta_gradient;

    % sanity check
    % if isnan(gradientUpdate)
    %    fprintf('gradient is NAN. Quit.\n');
    %    break;
    % end
   
    % gradient descent
    % theta = theta + gradientUpdate';

    % projection of all parameters to the feasible set
    % the feasible set of parameters in this case is greater than 0.1
    % if any(theta < 0.1)
    %    neg_indicator = (theta < 0.1);
    %    pos_indicator = ~neg_indicator;
    %    theta_min = 0.1*ones(4,1);
    %    theta = neg_indicator.*theta_min + pos_indicator.*theta;
    % end

    % store the parameters
    % param_hist = [param_hist theta];

    % terminate if the total number of iterations is more than maxIterations
    % if itr >= maxIterations
    %    break;
    % end
% end


%% plot trajectory

%% Debug session
% check_dx_dtheta = sum(isnan(dx_dtheta),'all');
% check_du_dtheta = sum(isnan(du_dtheta),'all');
