function lik = likBern(Bern,zPos,lik_matrix,model,filter)
%Calculate association likelihood for a bernoulli component with Gaussian
%process representation
num_meas = size(zPos,2);
if model.occlusion
    useNI = Bern.gammaCorrectionFactor == 1; %Might be better to put the whole struct into the filter.
    likGP = filter.calculateLogLikelihoodNegative(Bern.GPState(end).state,Bern.GPState(end).cov,zPos,lik_matrix,useNI);
    %likGP = filter.measLogLikelihoodNegative(Bern.GPState(end).state,Bern.GPState(end).cov,zPos,useNI);
    correctionFactor = Bern.gammaCorrectionFactor;
    [~,~,likGamma] = updateGamma(Bern.GPState(end).alpha,Bern.GPState(end).beta,num_meas,correctionFactor);
    lik = likGP + likGamma + log(Bern.r) + log(Bern.Pd) + log(Bern.w_death(end));
else
    likGP = filter.calculateLogLikelihood(Bern.GPState(end).state,Bern.GPState(end).cov,zPos,lik_matrix);
    %likGP2 = filter.measLogLikelihood(Bern.GPState(end).state,Bern.GPState(end).cov,zPos);
    [~,~,likGamma] = updateGamma(Bern.GPState(end).alpha,Bern.GPState(end).beta,num_meas,1);
    lik = likGP + likGamma + log(Bern.r) + log(model.Pd) + log(Bern.w_death(end));
end
end
