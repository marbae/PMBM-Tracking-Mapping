function [Bern,lik] = detectionPPP(zPos,PPP,indices,model,filter)
 
wp = PPP.w(indices);
GPStates = PPP.GPState(indices);
if model.occlusion
    Pd = PPP.Pd(indices);
    gammaCorrectionFactor = PPP.gammaCorrectionFactor(indices);
end

nComponents = length(wp);
num_meas = size(zPos,2);
for i = 1:nComponents
    if model.occlusion
        filter.useNI = gammaCorrectionFactor(i) == 1; %Might be better to put the whole struct into the filter.;
        [GPStates(i).state,GPStates(i).cov,lik_GP(i)] = filter.update(GPStates(i).state,GPStates(i).cov,zPos,model);
        [GPStates(i).alpha,GPStates(i).beta,likGamma(i)] = updateGamma(GPStates(i).alpha,GPStates(i).beta,num_meas,gammaCorrectionFactor(i));
    else
        filter.useNI = 0;
        [GPStates(i).state,GPStates(i).cov,lik_GP(i)] = filter.update(GPStates(i).state,GPStates(i).cov,zPos,model);
        [GPStates(i).alpha,GPStates(i).beta,likGamma(i)] = updateGamma(GPStates(i).alpha,GPStates(i).beta,num_meas,1);
    end
end
if model.occlusion
    w_c = lik_GP' + likGamma' + wp + log(Pd);
else
    w_c = lik_GP' + likGamma' + wp + log(model.Pd);
end
[w_hat,Bern.GPState] = GGP_merge(w_c,GPStates,model);
%Note: w_hat is not log-valued so this is done here
if num_meas > 1
    lik = log(w_hat);
    Bern.r = 1;
else
    Bern.r = w_hat/(w_hat+model.lambda_fa);
    lik = log(w_hat+model.lambda_fa);
end

% lik = log(w_hat + model.lambda_fa^num_meas);


end

