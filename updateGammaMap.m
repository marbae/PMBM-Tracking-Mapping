function [updAlpha,updBeta,lik] = updateGammaMap(alpha,beta,num_meas,correctionFactor,model)
%UPDATEGAMMA Update step of a bayesian estimation of the gamma distribution
%   Detailed explanation goes here

[alpha,beta] = predictGamma(alpha,beta,model.eta);

alpha = correctionFactor*alpha; %Corrects for known partial occlusion
% beta = correctionFactor*beta;
%num_meas = num_meas/correctionFactor;
if(correctionFactor == 0)
    alpha = alpha+1;
end
updAlpha = alpha+num_meas;
updBeta = beta+1;
lik = gammaln(updAlpha)-gammaln(alpha)+alpha*log(beta)-updAlpha*log(updBeta)-gammaln(num_meas+1);
if correctionFactor>0
    updAlpha = updAlpha/correctionFactor; %Recorrects to what would be the normal expected number of measurements 
 %   updBeta = updBeta/correctionFactor;
end

