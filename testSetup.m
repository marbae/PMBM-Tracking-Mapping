function [model,zPos,zPosFull,groundTruth] = testSetup(testFile)


model.sigma_v = 0.2;  %standard deviation of motion noise
model.sigma_yaw = 0.1;
model.sigma_r = 0.5;  %standard deviation of measurement noise

w_e_gamma = 20;
model.eta = 1/(1-1/w_e_gamma);

%range of the surveillance area
range_c = [-1 1;-1 1]*50;
%Poisson false alarm (clutter) rate
lambda_c = 60; %Number of clutter measurements per time step? m_k
model.sensorRes = deg2rad(1);
%Poisson clutter intensity
model.lambda_fa = lambda_c/prod(range_c(:,2)-range_c(:,1));
model.Nkin = 6;

%re-construct the data structure
testFile = "C:\Users\martibae\Documents\GP_PMBM\Data\" + testFile + ".mat"; %Location of mat files
load(testFile,'zPos','zPosFull','groundTruth')

dt = 0.1;
% generate tracks (ground truth) %Can loop this for several generated
% tracks, something similar for measurements.
L = 3.8;
W = 1.6;
D = 1.3;
S = 1.6;
N_theta_test = 9;

s = rng(1837);

model.K = length(zPos);
model.Ts = dt;   %sampling interval
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
model.N_theta_test=N_theta_test;
%target existence probability
model.Ps = 0.99;
%target detection probability
model.Pd = 0.99;
if(~isa(groundTruth,"cell"))
    groundTruthcell = cell(1,1);
    groundTruth.tbirth = 200;
    groundTruth.tdeath = 1220;
    groundTruthcell{1} = groundTruth;
    groundTruth = groundTruthcell;
else
    for i = 1:length(groundTruth)
        groundTruth{i}.tbirth = 950;
        groundTruth{i}.tdeath = 1575;
    end
end


% target initial state
model.vesselLength=4; %7 or 4
model.vesselWidthRatio = 2;
model.vesselBowRatio = 2;
model.initAlpha = 1000; %500 or 1000
model.initBeta = 100;
model.birthcov = blkdiag(diag([20,20,pi/2,3,3,pi/4]),model.gp.Ktt);
model.v0=2;

%Generate generic birth density
range = 40; %40 or 65
nBirths = 36;
model.birth = GenerateBirthPPP1(model,range,nBirths);

% if(multitarget)
%     startIndex1 = 500;
%     startIndex2 = 525;
%     nbirths = 2;
%     xstart = zeros(6,nbirths); %motion model dimensionality
%     diffx1 = ground1(1,startIndex1+10)-ground1(1,startIndex1);
%     diffy1 = ground1(2,startIndex1+10)-ground1(2,startIndex1);
%     heading1 = atan2(diffy1,diffx1);
%     xstart(:,1) =  [ground1(1,startIndex1); ground1(2,startIndex1); heading1; diffx1; diffy1; 0];
%     diffx2 = ground2(1,startIndex2+10)-ground2(1,startIndex2);
%     diffy2 = ground2(2,startIndex2+10)-ground2(2,startIndex2);
%     heading2 = atan2(diffy2,diffx2);
%     xstart(:,2) =  [ground2(1,startIndex2); ground2(2,startIndex2); heading2; diffx2; diffy2; 0];
% else
%     nbirths = 1;
%     xstart = zeros(6,nbirths); %motion model dimensionality
%     startState = targetcorr{1,770};
%     xstart(:,1) =  [startState.N; startState.E; deg2rad(startState.heading); startState.velN; startState.velE; deg2rad(startState.angVel)];
% end

% target initial state
% nbirths = length(groundTruth);
% xstart = zeros(model.Nkin,nbirths);
% startIndices = [500 525];
% for i = 1:nbirths
%     if(nbirths>1)
%         pos = groundTruth{i}.xKin(:,startIndices(i));
%         extend = groundTruth{i}.extendRadii;
%     else
%         pos = groundTruth.xKin(:,200);
%         extend = groundTruth.extendRadii;
%     end
%     xstart(:,i) = pos;
%     xstart(6,i) = deg2rad(pos(6));
%     xExtendEst0(:,i)=extend(1:360/N_theta_test:end);
% end
%xstart(:,1) = [-40; 12; 0; 0; 0; 0];
% xstart(:,2) = [27; -22; 3*pi/4; 0; 0; 0];
% xstart(:,1) = [-80 14 7*pi/4 0 0 0];
% xstart(:,3) = [75 75 0 0];
% xstart(:,4) = [75 -75 0 0];

%Birth model
% model.birth.w = 1/nbirths*ones(nbirths,1);
% for i = 1:nbirths
%     model.birth.GPState(i,1).state = [xstart(:,i); xExtendEst0(:,i)];
%     model.birth.GPState(i,1).cov = blkdiag(diag([20,20,pi,1,1,pi/20]),model.gp.Ktt);
%     model.birth.GPState(i,1).alpha = 10;
%     model.birth.GPState(i,1).beta = 2;
% end

% Gating parameters
Pg = 0.95;
model.gamma= chi2inv(Pg,2);
model.Qd = 1 - model.Pd*Pg;

% Thresholds
model.threshold_r = 1e-2;   %existence probability of Bernoulli component
model.threshold_u = 1e-2;   %weight of mixture component in PPP
model.threshold_w = 1e-2;   %weight of global hypothesis (multi-Bernoulli)
model.threshold_s = 1e-2;   %weight of the trajectory is still alive
model.recycle = 1e-1;       %recycling threshold
model.merge = 6;            %merge threshold used to merge similar GPs
model.M = 100;              %cap of number of MBM components in PMBM
model.num_iterations = 3;   %controls the number of iterations used in SO
model.max_repetition = 1;   %controls the number of iterations used in SO
model.maxIterGP = 50;       %Maximum GN iterations for Update step
model.GPMerge = false;      %Specific GP-merge strategy of Michaelis
model.GPvirtual = true;     %Whether to use virtual measurements
model.occlusion = true;     %Whether to modify model parameters to compensate for occlusion
model.positionCorrection = true;
model.headingCorrection = true;     %Whether to attempt to correct heading based on available measurements.
model.observatilityCriteria = false;


%extract target state from Bernoulli components with existence probability
%larger than this threshold 
model.exist_r = 0.5;        
