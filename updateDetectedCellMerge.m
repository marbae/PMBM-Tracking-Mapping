function [cell,lik] = updateDetectedCellMerge(cell,detectionMatrix,lik,model)
    nHypotheses = size(detectionMatrix,1);
    nDetections = sum(detectionMatrix,2);
    lik = zeros(1,nHypotheses);
    updAlpha = zeros(1,nHypotheses);
    updBeta = zeros(1,nHypotheses);
    [unCell,likFalseDetection] = updateUndetectedCell(cell,model);
    wFalse = 0;
    wTrue = 0;
    for i = 1:nHypotheses
        if(nDetections(i) > 0)
            [updAlpha(i),updBeta(i),likGamma] = updateGammaMap(cell.alpha,cell.beta,nDetections(i),cell.corrFactor,model);
            lik(i) = cell.occProb*cell.Pd*exp(likGamma)*nDetections(i)/model.delta^2;
            wTrue = lik(i)+wTrue;
        else
            updAlpha(i) = unCell.alpha;
            updBeta(i) = unCell.beta;
            lik(i) = likFalseDetection;
            wFalse = lik(i)+wFalse;
        end
    end
    wFalse = wFalse/sum(lik);
    wTrue = wTrue/sum(lik);
    lik = lik/sum(lik);
    
    [cell.alpha,cell.beta] = gammaMerge(updAlpha,updBeta,lik);
    cell.occProb = 0.9999*(wTrue+wFalse*unCell.occProb);
    lik = sum(lik);
end