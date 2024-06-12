function [w_hat,GPState] = GGP_merge(w_c,GPStates,model)
%GGP_MERGE Merging (mixture reduction) of a gamma GP distribution
%   Detailed explanation goes here
nComponents = size(GPStates,1);

for i = 1:nComponents
    means(:,i) = GPStates(i).state;
    covs(:,:,i) = GPStates(i).cov;
    alphas(i) = GPStates(i).alpha;
    betas(i) = GPStates(i).beta;
end
w = exp(w_c);
if(any(w==0))
    w(:) = 1e-100; %If the likelihood is so low as to 0 exp
end
wsum = sum(w);
n_x = size(means,1); 
if model.GPMerge %(Michaelis 2019 merge strategy)
    newState = sum(means(1:model.Nkin,:).*repmat(w(:)',model.Nkin,1),2)/wsum;
    newState(3) = angleMean(means(3,:),w);
    newRadii = Extentmerge(nComponents, means, model, newState);
    newExtent = sum(newRadii.*repmat(w(:)',n_x-model.Nkin,1),2)/wsum;
    newMean = [newState; newExtent];
else
    n_x = size(means,1); 
    newMean = sum(means.*repmat(w(:)',n_x,1),2)/wsum;
    newMean(3) = angleMean(means(3,:),w);
end
newCov = zeros(n_x,n_x);
for i = 1:nComponents
    dMean = means(:,i)-newMean;
    dMean(3) = angleDistance(means(3,i),newMean(3));
    newCov = newCov + w(i)*(covs(:,:,i)+(dMean)*(dMean)');
end
newCov = newCov/wsum;
GPState.state = newMean;
GPState.cov = newCov;

%Gamma merge (Granström and Orguner (2012))
[newAlpha,newBeta] = gammaMerge(alphas,betas,w);
GPState.alpha = newAlpha;
GPState.beta = newBeta;
w_hat = wsum;
end

