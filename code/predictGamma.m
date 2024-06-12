function [predAlpha,predBeta] = predictGamma(alpha,beta,eta)
%PREDICTGAMMA Prediction step of a bayesian estimation of a gamma
%distribution
%   Detailed explanation goes here
predAlpha = alpha/eta;
predBeta = beta/eta;
end

