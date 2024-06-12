function [w_hat,GPState] = Gaussian_merge(w_c,GPStates)

nComponents = size(GPStates,2);

for i = 1:nComponents
    means(:,i) = GPStates(i).state;
    covs(:,:,i) = GPStates(i).cov;
end

w = exp(w_c);
wsum = sum(w);
n_x = size(means,1); 
newMean = sum(means.*repmat(w(:)',n_x,1),2)/wsum;
newCov = zeros(n_x,n_x);
for i = 1:nComponents
    newCov = newCov + w(i)*(covs(:,:,i)+(means(:,i)-newMean)*(means(:,i)-newMean)');
end
newCov = newCov/wsum;
GPState.state = newMean;
GPState.cov = newCov;

w_hat = wsum;