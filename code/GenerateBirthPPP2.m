function birth = GenerateBirthPPP2(model,range,amount)

N_theta_test = model.N_theta_test;
dangle = deg2rad(360/N_theta_test);
L = model.vesselLength;
W = L/model.vesselWidthRatio;
D = L/model.vesselBowRatio;
S = W;
points = [range 0;0 range; -range 0; 0 -range];
shipExtend = ShipExtend(L, W, D, S, 'boxParabolicBow', dangle);
angles = model.gp.theta_test;
nbirths = 36;
%angles = deg2rad(linspace(0,360,amount));
xstart = zeros(model.Nkin,nbirths);
v0=model.v0;
ind = 1;
for i = 1:4
    for j = 1:9
        state = zeros(6,1);
        state(1:2) = points(i,:)';
        state(3) = angles(j);
        [state(4),state(5)] = pol2cart(atan2(points(i,2),points(i,1)),v0);
        state(6) = 0;
        xstart(:,ind) = state;
        xExtendEst0(:,ind)=shipExtend.hull_r;
        ind = ind +1;
    end
end

%Birth model
birth.w = 1/nbirths*ones(nbirths,1);
for i = 1:nbirths
    birth.GPState(i,1).state = [xstart(:,i); xExtendEst0(:,i)];
    birth.GPState(i,1).cov = model.birthcov;
    birth.GPState(i,1).alpha = model.initAlpha;
    birth.GPState(i,1).beta = model.initBeta;
end
