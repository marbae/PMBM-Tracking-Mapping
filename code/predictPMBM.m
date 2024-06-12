function [PPP,MBM] = predictPMBM(PPP,MBM,model,birth,filter)

% Predict existing PPP
PPP.w = PPP.w + log(model.Ps);
PPPlen = length(PPP.GPState);
[state,cov] = arrayfun(@(x) filter.predict(x.state,x.cov,model.Ts), PPP.GPState,'UniformOutput',false);
[alpha,beta] = arrayfun(@(x) predictGamma(x.alpha,x.beta,model.eta), PPP.GPState,'UniformOutput',false);
for i = 1:PPPlen
    PPP.GPState(i).state = state{i};
    PPP.GPState(i).cov = cov{i};
    PPP.GPState(i).alpha = alpha{i};
    PPP.GPState(i).beta = beta{i};
end
% Incorporate PPP birth
PPP.w = [PPP.w;log(birth.w)];
PPP.GPState = [PPP.GPState;birth.GPState];

% Predict MBM
n_track = length(MBM.track);
for i = 1:n_track
    nh = length(MBM.track{i});
    for h = 1:nh
        if MBM.track{i}(h).Bern.w_death(end) >= model.threshold_s
            [newGPState.state,newGPState.cov] = filter.predict(MBM.track{i}(h).Bern.GPState(end).state,MBM.track{i}(h).Bern.GPState(end).cov,model.Ts);
            [newGPState.alpha,newGPState.beta] = predictGamma(MBM.track{i}(h).Bern.GPState(end).alpha,MBM.track{i}(h).Bern.GPState(end).beta,model.eta);
            MBM.track{i}(h).Bern.GPState = [MBM.track{i}(h).Bern.GPState;newGPState];
            MBM.track{i}(h).Bern.t_death = [MBM.track{i}(h).Bern.t_death MBM.track{i}(h).Bern.t_death(end)+1];
            MBM.track{i}(h).Bern.w_death = [MBM.track{i}(h).Bern.w_death(1:end-1) MBM.track{i}(h).Bern.w_death(end)*(1-model.Ps) MBM.track{i}(h).Bern.w_death(end)*model.Ps];
            %MBM.track{i}(h).Bern.r = MBM.track{i}(h).Bern.r*model.Ps;
        else
            MBM.track{i}(h).Bern.w_death(end) = 0;
        end
    end
end

end

