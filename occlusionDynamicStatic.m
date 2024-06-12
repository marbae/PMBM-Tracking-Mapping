function [localMap] = occlusionDynamicStatic(MBM,time,filter,model,localMap,egoState,cellsToUpdate)
    %Based on MBM, correct detection probability and gamma parameters based
    %on one target occluding another
    nHypotheses = length(MBM.w);
    N = egoState.N;
    E = egoState.E;
    for nh = 1:nHypotheses
       track_indices = find(MBM.table(nh,:)>0);
       nTar = length(track_indices);
       angles = zeros(nTar,2);
       varAngles = zeros(nTar,2);
       radii = zeros(nTar,1);
       bearingMax = zeros(nTar,1);
       bearingMin = zeros(nTar,1);
       existProb = zeros(nTar,1);
       for iTar = 1:nTar
           currTrack = MBM.track{track_indices(iTar)}(MBM.table(nh,track_indices(iTar)));
           if currTrack.Bern.t_death(end) == time && currTrack.Bern.w_death(end) > 0
               state = currTrack.Bern.GPState(end).state;
               cov = currTrack.Bern.GPState(end).cov;
               state(1) = state(1)-N;
               state(2) = state(2)-E;
               [angles(iTar,:),varAngles(iTar,:)] = filter.getVisibilityAngle(state,cov);
               bearingMax(iTar) = angles(iTar,2);
               bearingMin(iTar) = angles(iTar,1);
               radii(iTar,:) = sqrt(state(1).^2+state(2).^2);
               existProb(iTar) = currTrack.Bern.w_death(end);
           end
       end
    end
    if(nHypotheses>0)
        for i = 1:size(cellsToUpdate,1)
            ind = cellsToUpdate(i,:);
            currCell = localMap{ind(1),ind(2)}{1};
            r = sqrt((N-currCell.N)^2+(E-currCell.E)^2);
            deltaYmax = currCell.E-E+0.5*model.delta;
            deltaXmax = currCell.N-N+0.5*model.delta;
            deltaYmin = currCell.E-E-0.5*model.delta;
            deltaXmin = currCell.N-N-0.5*model.delta;
            [bMax,bMin] = getCellBearings(deltaYmax,deltaYmin,deltaXmin,deltaXmax);
            rCond = r>radii';
            bMaxCond = bMax<bearingMax;
            bMinCond = bMin>bearingMin;
            occlCond = rCond & bMaxCond & bMinCond;
            iOccl = find(occlCond);
            if(~isempty(iOccl))
                pNonOccl = 1;
                for k = 1:size(iOccl,2)
                    exist = existProb(iOccl(k));
                    pNonOccl = pNonOccl*(1-exist);
                end
                currCell.Pd = currCell.Pd*pNonOccl;
            end
            localMap{ind(1),ind(2)}{1} = currCell;
        end
    end
end


