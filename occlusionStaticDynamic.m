function [PPP,MBM] = occlusionStaticDynamic(MBM,PPP,egoState,localMap,time,filter,model)
    nHypotheses = length(MBM.w);
    mapSize = size(localMap,2);
    refCell = localMap{round(mapSize/2),round(mapSize/2)}{1};
    refPoint = [refCell.N; refCell.E];
    egoN = egoState.N;
    egoE = egoState.E;
    for nh = 1:nHypotheses
       track_indices = find(MBM.table(nh,:)>0);
       nTar = length(track_indices);
       angles = zeros(nTar,2);
       varAngles = zeros(nTar,2);
       radii = zeros(nTar,1);
       existProb = zeros(nTar,1);
       for iTar = 1:nTar
           currTrack = MBM.track{track_indices(iTar)}(MBM.table(nh,track_indices(iTar)));
           if currTrack.Bern.t_death(end) == time && currTrack.Bern.w_death(end) > 0
               state = currTrack.Bern.GPState(end).state;
               cov = currTrack.Bern.GPState(end).cov;
               %[angles(iTar,:),varAngles(iTar,:)] = filter.getVisibilityAngle(state,cov,egoState);
               pos = [state(1) state(2)];
               mapPos = pos-refPoint';
               mapPos = round(mapPos/model.delta);
               r = sqrt(pos(1)^2+pos(2)^2);
                indexN = round(mapSize/2)+mapPos(1);
                indexE = round(mapSize/2)+mapPos(2);
                currTrack.Bern.Pd = localMap{indexN,indexE}{1}.Pd;
                if(currTrack.Bern.Pd == 0)
                    currTrack.Bern.Pd = 0.05;
                end
               MBM.track{track_indices(iTar)}(MBM.table(nh,track_indices(iTar))) = currTrack;
           end
       end
    end
    nPPP = length(PPP.w);
    PPP.Pd = model.Pd*ones(nPPP,1);
    nPPP = length(PPP.w);
    iToRemove = [];
    for iTar = 1:nPPP
        currComponent = PPP.GPState(iTar);
        state = currComponent.state;
        cov = currComponent.cov;
        pos = [state(1) state(2)];
        mapPos = pos'-refPoint;
        mapPos = round(mapPos/model.delta);
        r = sqrt(pos(1)^2+pos(2)^2);
        indexN = round(mapSize/2)+mapPos(1);
        indexE = round(mapSize/2)+mapPos(2);
        PPP.Pd(iTar) = localMap{indexN,indexE}{1}.Pd;
        if(PPP.Pd(iTar) < 0.05)
            iToRemove = [iToRemove iTar];
            PPP.Pd(iTar) = 0.05; %Consider using this to remove the PPP components that are irrelevant
        end
    end
    PPP.Pd(iToRemove) = [];
    PPP.GPState(iToRemove) = [];
    PPP.w(iToRemove) = [];
end

