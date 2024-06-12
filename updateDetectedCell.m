function [updAlpha,updBeta,updOccprob,lik] = updateDetectedCell(M,ind,nDetections,model)
    [updAlpha,updBeta,likGamma] = updateGammaMap(M.alpha(ind),M.beta(ind),nDetections,M.corrFactor(ind),model);
    likDetection = M.occProb(ind)*M.Pd(ind)*exp(likGamma)*nDetections/model.delta^2;
    [unAlpha,unBeta,unOccProb,likFalseDetection] = updateUndetectedCell(M,ind,model);
    likFalseDetection = likFalseDetection*model.lambda_fa^nDetections;
    wTrue = likDetection/(likDetection+likFalseDetection);
    wFalse = likFalseDetection/(likDetection+likFalseDetection);
    [updAlpha,updBeta] = gammaMerge([updAlpha;unAlpha],[updBeta;unBeta],[wTrue;wFalse]);
    updOccprob = 0.9999*(wTrue+wFalse*unOccProb);
    lik = 1;
end