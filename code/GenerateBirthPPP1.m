function birth = GenerateBirthPPP1(model,egopos)

egoN = egopos.N;
egoE = egopos.E;

N_theta_test = model.N_theta_test;
dangle = deg2rad(360/N_theta_test);
L = model.vesselLength;
W = L/model.vesselWidthRatio;
D = L/model.vesselBowRatio;
S = W;
shipExtend = ShipExtend(L, W, D, S, 'boxParabolicBow', dangle);
nbirths = model.nBirths;
angles = deg2rad(linspace(0,360,nbirths));
xstart = zeros(model.Nkin,nbirths);
v0=model.v0;
for i = 1:nbirths
    state = zeros(6,1);
    [state(1),state(2)] = pol2cart(angles(i),model.range);
    state(1) = state(1)+egoN;
    state(2) = state(2)+egoE;
    state(3) = mod(angles(i)-pi,2*pi);
    [state(4),state(5)] = pol2cart(state(3),v0);
    state(6) = 0;
    xstart(:,i) = state;
    xExtendEst0(:,i)=shipExtend.hull_r;
end

%Birth model
birth.w = 1/nbirths*ones(nbirths,1);
for i = 1:nbirths
    birth.GPState(i,1).state = [xstart(:,i); xExtendEst0(:,i)];
    birth.GPState(i,1).cov = model.birthcov;
    birth.GPState(i,1).alpha = model.initAlpha;
    birth.GPState(i,1).beta = model.initBeta;
end
