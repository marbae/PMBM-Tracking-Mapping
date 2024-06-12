function [model,fullZ,groundTruth,egoState,environment] = simulationSetup(simulationFile)

%simulationFile = "C:\Users\martibae\Documents\GP_PMBM\Simulations\" + simulationFile + ".mat";
simulationFile = simulationFile + ".mat";
%Groundtruth is a kinematic state and a star-convex shape, Scenario is
%measurements and probability of detection along with clutter rate.
load(simulationFile,'groundTruth','Scenario','egoState','environment')

zPos = Scenario.Z{1};
fullZ = Scenario.Z;
K = length(zPos);
model.K = K;

% Effective window length for the gamma prediction
w_e_gamma = 20;

model.sigma_v = 0.2;  %standard deviation of motion noise
model.sigma_yaw = 0.1;
model.sigma_r = 0.3;  %standard deviation of measurement noise
model.eta = 1/(1-1/w_e_gamma);
model.Ts = 1;   %sampling interval

%range of the surveillance area
range_c = [-1 1;-1 1]*100;
%Poisson false alarm (clutter) rate
lambda_c = 60;
%Poisson clutter intensity
%model.lambda_fa = lambda_c/prod(range_c(:,2)-range_c(:,1));
model.lambda_fa = Scenario.clutter_intensity;
model.sensorRes = deg2rad(0.5);

%target existence probability
model.Ps = 0.95;
%target detection probability
model.Pd = 0.99;
%re-construct the data structure

% generate tracks (ground truth) %Can loop this for several generated
% tracks, something similar for measurements.
s = rng(1837);
N_theta_test = 9;
model.Nkin=6;

%Generate measurements
sigma_f=0.5; %Kernel magnitude 
sigma_r=0.5; %Baseline radius uncertainty
sigma_n=1e-2; %Nugget
l=pi/4; %Scale (how far away are angles correlated)
hyperparameters = [sigma_f sigma_r sigma_n l]';
dtheta_plot = deg2rad(10);
forgetting_factor = 1e-2;
symmetric = true;
model.gp = GP(hyperparameters, symmetric, N_theta_test, dtheta_plot, forgetting_factor);
model.Nx = N_theta_test+model.Nkin;
model.N_theta_test = N_theta_test;
% target initial state
nbirths = length(groundTruth);
xstart = zeros(model.Nkin,nbirths);
% for i = 1:nbirths
%     pos = groundTruth{i}.xKin(:,70);
%     extend = groundTruth{i}.extendRadii;
%     xstart(:,i) = pos;
%     xExtendEst0(:,i)=extend(1:360/N_theta_test:end);
% end

%Birth model for generic birth
model.vesselLength=6;
model.vesselWidthRatio = 2;
model.vesselBowRatio = 3;
model.initAlpha = 500;
model.initBeta = 100;
model.birthcov = blkdiag(diag([20,20,pi/2,3,3,pi/4]),model.gp.Ktt);
model.v0=2;

%Generate generic birth density
model.range = 105;
model.nBirths = 36;
%model.birth.w = model.birth.w*0.10;

% Gating parameters
Pg = 0.99;
model.gamma= chi2inv(Pg,2);
model.Qd = 1 - model.Pd*Pg;

%Mapping parameters
model.maxRange = model.range-5;
model.delta=1;
model.initProb = 0.5;

% Thresholds
model.threshold_r = 1e-2;   %existence probability of Bernoulli component (not used)
model.threshold_u = 1e-2;   %weight of mixture component in PPP
model.threshold_w = 1e-2;   %weight of global hypothesis (multi-Bernoulli)
model.threshold_s = 1e-1;   %weight of the trajectory is still alive
model.recycle = 1e-1;       %recycling threshold (not used)
model.merge = 4;            %merge threshold used to merge similar GGIWs
model.M = 100;              %cap of number of MBM components in PMBM
model.num_iterations = 3;   %controls the number of iterations used in SO
model.max_repetition = 1;   %controls the number of iterations used in SO
model.maxIterGP = 50;       %Maximum GN iterations for Update step

%Toggles for extra features (occlusion is mandatory for mapping)
model.GPMerge = false;      %Whether to use the specific GP-merge strategy
model.GPvirtual = false;     %Whether to use virtual measurements
model.occlusion = true;     %Whether to modify model parameters to compensate for occlusion
model.positionCorrection = true;
model.headingCorrection = true;     %Whether to attempt to correct heading based on available measurements.
model.observatilityCriteria = false;

%extract target state from Bernoulli components with existence probability
%larger than this threshold 
model.exist_r = 0.5;        
