function [Bern,lik] = misdetectionBern(Bern,model)
if model.occlusion
    Pd = Bern.Pd;
else
    Pd = model.Pd;
end
Qd = 1-Pd;
alpha=Bern.GPState(end).alpha;
beta=Bern.GPState(end).beta;
zeroDetection = Bern.w_death(end)*Pd*(beta/(beta+1))^alpha; %Likelihood of cardinality 0 from the poisson process
misDetection = 1-(1-Qd)*Bern.w_death(end); %Likelihood of missing detection from the detection model

qD = misDetection + zeroDetection; %total likelihood of not detecting a target

w1 = misDetection/qD;
w2 = zeroDetection/qD;
lik = 1 - Bern.r + Bern.r*qD;

[Bern.GPState(end).alpha,Bern.GPState(end).beta] = gammaMerge([alpha;alpha],[beta;beta+1],[w1;w2]);

Bern.r = Bern.r*qD/lik;
lik = log(lik);
Bern.w_death = [Bern.w_death(1:end-1) Bern.w_death(end)*(qD)]/(1-Bern.w_death(end)*(1-qD)); %Termination component
end
