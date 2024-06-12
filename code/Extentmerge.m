function newRadii = Extentmerge(nComponents, means, model, newState)
%Calculates new radii with the new merged state for the components in the
%mixture by performing GP regression
    for i = 1:nComponents
        oldMean = means(:,i);
        [localPoints(1,:),localPoints(2,:)] = pol2cart(model.gp.theta_test',oldMean(model.Nkin+1:end));
        globalPoints = oldMean(1:2)+rotate(localPoints,oldMean(3));
        newLocalPoints = globalPoints - newState(1:2);
        newLocalPoints = rotate(newLocalPoints,-newState(3));
        [testAngles,testRadii] = cart2pol(newLocalPoints(1,:),newLocalPoints(2,:));
        if(any(testRadii<0))
            Ktt = model.gp.covarianceMatrix(testAngles,testAngles);
            Kzt = model.gp.covarianceMatrix(model.gp.theta_test,testAngles);
            Hzt = Kzt/Ktt;
            HRadii = Hzt;
            newRadii(:,i) = HRadii*testRadii';
        else
            newRadii(:,i) = testRadii';
        end
    end
end