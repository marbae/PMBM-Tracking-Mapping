clear all
close all
load("egoSim_1Tracks.mat")
zPos = Scenario.Z{1};
model.maxRange = 100;
model.delta=1;
model.initProb = 0.5;
model.lambda_fa = 1e-4;
model.Pd = 0.95;
w_e_gamma = 20;
model.eta = 1/(1-1/w_e_gamma);
fullMap = model.initProb*ones(round(1000/model.delta));
measMap = 80*ones(round(1000/model.delta));
localMapSize = round(model.maxRange*3/model.delta);
localMap.alpha = 800*model.delta*ones(localMapSize);
localMap.beta = 10*ones(localMapSize);
localMap.occProb = model.initProb*ones(localMapSize);
localMap.corrFactor = 1*ones(localMapSize);
localMap.Pd = model.Pd*ones(localMapSize);
refX = round(egoState{1}.N);
refY = round(egoState{1}.E);
tic
for i = 1:localMapSize
    for j = 1:localMapSize
        localMap.N(i,j) = round(refX+i*model.delta)-round(localMapSize/2)*model.delta;
        localMap.E(i,j) = round(refY+j*model.delta)-round(localMapSize/2)*model.delta;
    end
end

%Main loop
occupiedCells = zeros(localMapSize);
times = zeros(4,length(egoState));
mapVis.allOccMap = zeros(localMapSize,localMapSize,length(egoState));
mapVis.allN = mapVis.allOccMap;
mapVis.allE = mapVis.allN;
for t = 1:length(egoState)
    tic
    [localMap,cellsToUpdate] = calculateMapPd(egoState{t},localMap,occupiedCells,model); %Requires smart implementation to not query whole map many times
    times(1,t) = toc; 
    tic
    globalZ = convertLocalToGlobal(egoState{t},zPos{t});
    %Data association, only most simple for now.
    detectedCell = findAssocCell(globalZ,localMap,model); %Simple rounding based on the centerpoint of the map
    times(2,t) = toc;
    tic
    %Map Update
    alpha = localMap.alpha;
    beta = localMap.beta;
    occProb = localMap.occProb;
    for i = 1:size(cellsToUpdate,1)
        ind = cellsToUpdate(i,:);
        if(detectedCell(ind)>0)
            occupiedCells(ind) = 1; %Not sure about this
            [alpha(ind),beta(ind),occProb(ind),lik] = updateDetectedCell(localMap,ind,detectedCell(ind),model); %Simple rules
        else
            [alpha(ind),beta(ind),occProb(ind),lik] = updateUndetectedCell(localMap,ind,model); %Simple rules
        end
    end
    localMap.alpha = alpha;
    localMap.beta = beta;
    mapVis.allOccMap(:,:,t) = occProb;
    mapVis.allN(:,:,t) = localMap.N;
    mapVis.allE(:,:,t) = localMap.E;
    localMap.occProb = occProb;
    times(3,t) = toc;
    tic
    [localMap,fullMap,measMap,occupiedCells] = moveLocalMap(egoState{t},localMap,fullMap,measMap,occupiedCells,model); %Requires some thought to make it good
    times(4,t) = toc;
    disp(t)
%     imagesc(occProb)
%     pause(0.5)
end
%%
plotTimeLapseLocal(egoState,groundTruth,environment,zPos,mapVis,model,1,0,0.2,'Leuvensim1')




