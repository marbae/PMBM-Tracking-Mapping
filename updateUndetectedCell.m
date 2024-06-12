function [updAlpha,updBeta,updOccProb,lik] = updateUndetectedCell(M,ind,model)
    preAlpha = M.alpha(ind); preBeta = M.beta(ind); preR = M.occProb(ind); corrFactor = M.corrFactor(ind);
    [preAlpha,preBeta] = predictGamma(preAlpha,preBeta,model.eta);
    Qd = 1-M.Pd(ind);
    alpha=preAlpha*corrFactor;
    beta=preBeta;%*cell.corrFactor;
    zeroDetection = M.Pd(ind)*(beta/(beta+1))^alpha; %Likelihood of cardinality 0 from the poisson process
    misDetection = 1-(1-Qd); %Likelihood of missing detection from the detection model
    
    qD = misDetection + zeroDetection; %total likelihood of not detecting a target
    
    w1 = misDetection/qD;
    w2 = zeroDetection/qD;
    lik = 1 - M.occProb(ind) + M.occProb(ind)*qD;
    
    [updAlpha,updBeta] = gammaMerge([preAlpha;preAlpha],[preBeta;preBeta+corrFactor],[w1;w2]);
    
    updOccProb = preR*qD/lik;
end