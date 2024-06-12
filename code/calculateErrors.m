function metrics = calculateErrors(assignment,groundTruth,estimate,time,model)
    nEstimates = size(estimate.state,2);
    nTargets = length(groundTruth);
    metrics = [];
    if(nTargets<=nEstimates)
        iEnd = nTargets;
    else
        iEnd = nEstimates;
    end
    for i = 1:iEnd
        if(nTargets<=nEstimates)
            iEst=assignment(i);
            iGT = i;
        else
            iEst=i;
            iGT = assignment(i);
        end
        if(groundTruth{iGT}.tbirth<=time && groundTruth{iGT}.tdeath>=time)
            state=estimate.state(:,iEst);
            GPState.state = state;
            GPState.cov = estimate.cov(:,:,iEst);
            hull_n = getExtendCoordinates(GPState,model);
            metrics.iGT(i) = iGT;
            metrics.IOU(i)=IntersectOverUnion(hull_n',groundTruth{iGT}.extendNE(:,:,time)');
            metrics.RMSE(i) = norm(state(1:2)-groundTruth{iGT}.xKin(1:2,time));
            metrics.headingErr(i) = angleDistance(state(3),groundTruth{iGT}.xKin(3,time));%abs(asin(sin(ssa(state(3))-ssa(groundTruth{iGT}.xKin(3,time)))));
            metrics.RMSEvel(i) = norm(state(4:5)-groundTruth{iGT}.xKin(4:5,time));
            GTextent = groundTruth{iGT}.extendRadii(1:360/model.N_theta_test:end);
            GTkin = groundTruth{iGT}.xKin(:,time);
            GTstate = [GTkin; GTextent];
            delta = GTstate-state;
            delta(3) = angleDistance(state(3),GTstate(3));
            metrics.NEES(i) = delta'/(estimate.cov(:,:,iEst))*delta;
            vec = [1 2 4 5];
            metrics.NEESposvel(i) = delta(vec)'/(estimate.cov(vec,vec,iEst))*delta(vec);
            vec = [3 6];
            metrics.NEESheading(i) = delta(vec)'/(estimate.cov(vec,vec,iEst))*delta(vec);
            metrics.angvel(i) = delta(6);
        end
    end
end

