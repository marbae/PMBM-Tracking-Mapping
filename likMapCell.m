function lik = likMapCell(cell,zPos,model)
    nDetections = size(zPos,2);
    [~,~,likGamma] = updateGammaMap(cell.alpha,cell.beta,nDetections,cell.corrFactor,model);
    lik = log(cell.occProb)+log(cell.Pd)+likGamma+log(nDetections/model.delta^2);
end