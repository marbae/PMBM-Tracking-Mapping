function [ Error,metrics] = GOSPAmetric(estimate,groundTruth,time,c,model)

%
% Function that computes the Generalized Optimal Subpattern Assignment (OSPA) Metric
% for the random sets X and Y.
%
% Uses GP model for cost calculation
%
% Input X (and Y) is a struct with J estimates
% X.x -- Nx * J matrix with kinematic vectors
% X.X -- d * d * J tensor with random matrices

%Unpack input:

N_theta = model.N_theta_test;
p = 2;
Ny = length(groundTruth);
GTstate = zeros(model.Nx,Ny);
validTargets = false(1,Ny);
trueNy = 0;
for i = 1:Ny
%     GTextent = groundTruth{i}.extendRadii(1:360/N_theta:end);
      GTkin(1:2,i) = groundTruth{i}.xKin(1:2,time);
      GTkin(3:4,i) = groundTruth{i}.xKin(4:5,time);
      if(groundTruth{i}.tbirth<=time && groundTruth{i}.tdeath>=time)
          trueNy = trueNy+1;
          validTargets(i) = true;
     end
%       GTstate(:,i) = [GTkin; GTextent];
end
% Number of elements in X and Y.
Nx = size(estimate.state,2);
%Ny = trueNy;

% if there are no elements at all
if Nx==0 && Ny==0
    gospa = 0;
    LocationError = 0;
    MissedError = 0;
    FalseError = 0;
    Error = [gospa,LocationError,MissedError,FalseError];
    metrics = [];
    return
elseif Nx==0
    gospa = (0.5*(c^p)*trueNy)^(1/p);
    LocationError = 0;
    MissedError = trueNy;
    FalseError = 0;
    Error = [gospa,LocationError,MissedError,FalseError];
    metrics = [];
    return
elseif trueNy==0
    Nx = 0;
    gospa = (0.5*(c^p)*Nx)^(1/p);
    LocationError = 0;
    MissedError = 0;
    FalseError = Nx;
    Error = [gospa,LocationError,MissedError,FalseError];
    metrics = [];
    return
end

% Distance matrix
D = repmat(c,[Nx Ny]);

for ix = 1:Nx
    for iy = 1:Ny
        % Gaussian Wasserstein Distance for kinematics and extent
        %delta = GTstate(:,iy)-estimate.state(:,ix);
        %distance = delta'/(estimate.cov(:,:,ix))*delta;
        %distance = GTstate(:,iy)-estimate.state(:,ix)
        estKin(1:2) = estimate.state(1:2,ix);
        estKin(3:4) = estimate.state(4:5,ix);
        distance = norm(estKin(1:2)'-GTkin(1:2,iy),2);
        % Poisson rates
        %prd = abs(X.g(ix)-Y.g(iy));
        
        % Apply threshold c
        D(ix,iy) = distance;
    end
end

% Allocate memory
gospa = 0;
LocationError = 0;
absGamma = 0;


if Ny<=Nx
    % Compute assignment
    [Customer2Item,~] = auctionAlgorithm(-log(D)');
    
    % Iterate over true targets
    for iy = 1:Ny
        % Check if distance is small enough
        if D(Customer2Item(iy),iy) < c && validTargets(iy)
            % Location part of GOSPA
            gospa = gospa + D(Customer2Item(iy),iy)^p;
            % Location error
            LocationError = LocationError + D(Customer2Item(iy),iy);
            % Number of assignments
            absGamma = absGamma + 1;
        end
    end
    
else
    % Compute assignment
    [Customer2Item,~] = auctionAlgorithm(-log(D));
    
    % Iterate over estimates
    for ix = 1:Nx
        % Check if distance is small enough
        if D(ix,Customer2Item(ix)) < c && validTargets(Customer2Item(ix))
            % Location part of GOSPA
            gospa = gospa + D(ix,Customer2Item(ix))^p;
            % Location error
            LocationError = LocationError + D(ix,Customer2Item(ix));
            % Number of assignments
            absGamma = absGamma + 1;
        end
    end
    
end

gospa = (gospa + 0.5*(c^p)*(Nx+trueNy-2*absGamma))^(1/p);

% Missed detection error
MissedError = trueNy-absGamma;

% False alarm error
FalseError = Nx-absGamma;

Error = [gospa,LocationError,MissedError,FalseError];
metrics = calculateErrors(Customer2Item,groundTruth,estimate,time,model);


% if Nx>Ny % assume X contains fewer elements than Y
%     tmp = X;
%     X = Y;
%     Y = tmp;
%     Nx = size(X,2);
%     Ny = size(Y,2);
% end
% 
% % Construct distance matrix D
% D = repmat(c,[Ny Ny]);
% for i = 1:Nx
%     xi = X(1:2,i);
%     for j = 1:Ny
%         xj = Y(1:2,j);
%         % Euclidean norm
%         D(i,j) = min(c,norm(xi-xj));
%     end
% end
% 
% % Use the auction algorithm to compute the best assignments
% % auction maximizes the reward. Maximizing the negative distance is equal
% % to minimizing the positive distance.
% [~,Item2Customer] = auctionAlgorithm(-D);
% 
% %D = dg+dx+dX;
% 
% % allocate memory
% alphas = repmat(c,[1 Ny]);
% % compute the alphas
% for j = 1:Ny
%     if Item2Customer(j) <= Nx
%         alphas(j) = D(Item2Customer(j),j);
%     end
% end
% 
% % Compute the OSPA metric
% gospa   = (mean(alphas.^p))^(1/p);


