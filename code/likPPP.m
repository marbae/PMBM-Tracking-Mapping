function lik = likPPP(zPos,lik_matrix,PPP,indices,model,filter)
%Calculate association likelihood for the bernoulli components in a PPP with Gaussian
%process representation (used in Stochastic Optimization to find likelihood
%for use in the optimization)

wp = PPP.w(indices);
GPStates = PPP.GPState(indices);
lik_matrix = lik_matrix(:,indices);
if model.occlusion
    Pd = PPP.Pd(indices);
    gammaCorrectionFactor = PPP.gammaCorrectionFactor(indices);
end

nComponents = length(wp);
num_meas = size(zPos,2);
for i = 1:nComponents
    if model.occlusion
        useNI = gammaCorrectionFactor(i) == 1;
        lik_GP(i) = filter.calculateLogLikelihoodNegative(GPStates(i).state,GPStates(i).cov,zPos,lik_matrix(:,i),useNI);
        %lik_GP(i) = filter.measLogLikelihoodNegative(GPStates(i).state,GPStates(i).cov,zPos,useNI);
        [~,~,likGamma(i)] = updateGamma(GPStates(i).alpha,GPStates(i).beta,num_meas,gammaCorrectionFactor(i));
    else
        lik_GP(i) = filter.calculateLogLikelihood(GPStates(i).state,GPStates(i).cov,zPos,lik_matrix(:,i));
        %lik_GP(i) = filter.measLogLikelihood(GPStates(i).state,GPStates(i).cov,zPos);
        [~,~,likGamma(i)] = updateGamma(GPStates(i).alpha,GPStates(i).beta,num_meas,1);
    end

end
if model.occlusion
    w_c = lik_GP' + likGamma' + wp + log(Pd);
else
    w_c = lik_GP' + likGamma' + wp + log(model.Pd);
end
% % if(size(w_c,1)>1)
%     [w_hat,~] = GGP_merge(w_c,GPStates,model);
% % else
% %     w_hat = w_c;
% % end
w_hat = sum(exp(w_c));
num_meas = size(zPos,2);
if num_meas > 1
    lik = log(w_hat);
else
    lik = log(w_hat+model.lambda_fa);
end

%lik = log(w_hat + model.lambda_fa^num_meas);
end

