function [Bern,lik] = detectionBern(Bern,zPos,model,filter)
%Performs Kalman filter Update for one Bernoulli component
num_meas = size(zPos,2);
if model.occlusion
    filter.useNI = Bern.gammaCorrectionFactor == 1; %Might be better to put the whole struct into the filter.
    [Bern.GPState(end).state,Bern.GPState(end).cov,likGP] = filter.update(Bern.GPState(end).state,Bern.GPState(end).cov,zPos,model);
    correctionFactor = Bern.gammaCorrectionFactor;
    [Bern.GPState(end).alpha,Bern.GPState(end).beta,likGamma] = updateGamma(Bern.GPState(end).alpha,Bern.GPState(end).beta,num_meas,correctionFactor);
    lik = likGP + likGamma + log(Bern.r) + log(Bern.Pd) + log(Bern.w_death(end));
else
    [Bern.GPState(end).state,Bern.GPState(end).cov,likGP] = filter.update(Bern.GPState(end).state,Bern.GPState(end).cov,zPos,model);
    [Bern.GPState(end).alpha,Bern.GPState(end).beta,likGamma] = updateGamma(Bern.GPState(end).alpha,Bern.GPState(end).beta,num_meas,1);
    lik = likGP + likGamma + log(Bern.r) + log(model.Pd) + log(Bern.w_death(end));
end
Bern.r = 1;

Bern.t_death = Bern.t_death(end);
Bern.w_death = 1;

end
