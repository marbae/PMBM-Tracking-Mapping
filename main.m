clear;clc
dbstop if error
close all
warning('off','MATLAB:polyshape:repairedBySimplify') %IOU calculation throws this

dataSource = 'sim';
testFile = "2022-08-30-14-56-31";
simFile = "egoSim_0Tracks";
if(strcmp(dataSource,'sim'))
    [model,zPos,groundTruth,egoState,environment] = simulationSetup(simFile);
    numMC = 1;
    videofile = simFile+'_GP';
end
if(strcmp(dataSource,'test'))
    [model,zPos,zPosFull,groundTruth] = testSetup(testFile);
    numMC = 1;
    videofile = testFile+'_GP';
end

%Initialization of maps
fullMap = model.initProb*ones(round(1000/model.delta));
localMapSize = round(model.maxRange*3/model.delta);
localMap = cell(localMapSize);
mapCell.alpha = 400*model.delta;
mapCell.beta = 10;
mapCell.occProb = model.initProb;
mapCell.corrFactor = 1;
mapCell.Pd = 0.95;

%Parameters used in GOSPA metric
c = 10;
p = 2;

%Number of time steps
K = model.K;
GOSPA = zeros(K,4,numMC);
trajectoryEstimates = cell(numMC,1);
simulation_times = zeros(1,numMC);
estimates = cell(K,1);
metrics = cell(K,numMC);
%%
for t = 1:numMC
    simulation_time = zeros(K,1);
    if(strcmp(dataSource,'sim'))
        Z=zPos{t};
    else
        Z=zPos;
    end
    birth = GenerateBirthPPP1(model,egoState{1});
    % Initialisation
    PPP.w = log(birth.w); %General note: weights and likelihoods are stored in logarithmic form (Could be a source of error)
    PPP.GPState = birth.GPState;
    MBM.w = [];     % Global hypotheses weights
    MBM.track = {}; % Local hypotheses trees
    MBM.table = []; % Global hypotheses look-up table
    refX = round(egoState{1}.N);
    refY = round(egoState{1}.E);
    for i = 1:size(localMap,1)
        for j = 1:size(localMap,2)
            currCell = mapCell;
            currCell.N = round(refX+i*model.delta)-round(localMapSize/2)*model.delta;
            currCell.E = round(refY+j*model.delta)-round(localMapSize/2)*model.delta;
            localMap{i,j} = {currCell};
        end
    end
    occupiedCells = zeros(localMapSize);
    
    estimates = cell(K,1);
    trajectoryEstimates{t} = cell(K,1);
    cv3 = CV([model.sigma_v; model.sigma_v; model.sigma_yaw]);
    filter = ExtendedGP(cv3,model);
    times = zeros(7,K);
    for k = 1:K
        if(strcmp(dataSource,'sim'))
            meas = Z{k};
        else
            meas = Z{k}';
        end
        globalMeas = convertLocalToGlobal(egoState{k},meas);
        %Print info
        [t,k]
        %Update step
        [PPP,MBM,localMap,occupiedCells,times(:,k)] = updatePMBM(PPP,MBM,globalMeas,...
            k,model,filter,localMap,occupiedCells,egoState{k});
        %Extract estimates (both estimate of the current time and the
        %estimate of the full trajectory)
        [estimates{k},trajectoryEstimates{t}{k}] = estimator(MBM,model);
        length(trajectoryEstimates{t}{k})
        %Evaluate filtering performance using GOSPA
        [GOSPA(k,:,t),metrics{k,t}] = GOSPAmetric(estimates{k},groundTruth,k,c,model);
        %Prediction Step
        if k < K
            birth = GenerateBirthPPP1(model,egoState{k+1});
            [PPP,MBM] = predictPMBM(PPP,MBM,model,birth,filter);
        end
        %Construct a new local map when we move to far off center from the
        %old one
        [localMap,fullMap,occupiedCells] = moveLocalMap(egoState{k},localMap,fullMap,occupiedCells,model);
        times(7,k) = toc;
        simulation_time(k) = sum(times(:,k));
    end
    simulation_times(:,t) = sum(simulation_time);
end
%%
doPlot = 1;
doMovie = 1;
plotMC = 1;
if(doPlot)
    generatePlots(metrics,model,length(groundTruth),plotMC)
end
if(doMovie)
    plotTimeLapse(groundTruth, Z, trajectoryEstimates{numMC,1}, [-100 100 -100 100], model.Ts,1,0,videofile)
end
[GOSPAtable,metricsTable]=retrieveTable(real(GOSPA),metrics,length(groundTruth));
meanGOSPA = real(mean(GOSPAtable,2));
meanMetrics = mean(metricsTable,2);
f6 = figure('Name','Average Over MC runs');
figHandle = generateSubplot(metrics,real(GOSPA),'GP',[950 1550],f6,'-');
save(videofile,'GOSPA','metrics','model')
