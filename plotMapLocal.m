function [] = plotMapLocal(ego,occMap,N,E,model)
iPlot = find(occMap>0.55);
nPlot = size(iPlot,1);
for i = 1:nPlot
    ncoord = [N(iPlot(i))-model.delta/2 N(iPlot(i))-model.delta/2 N(iPlot(i))+model.delta/2 N(iPlot(i))+model.delta/2];
    ecoord = [E(iPlot(i))-model.delta/2 E(iPlot(i))+model.delta/2 E(iPlot(i))+model.delta/2 E(iPlot(i))-model.delta/2];
    ncoord = ncoord-ego.N;
    ecoord = ecoord-ego.E;
    theta = deg2rad(ego.heading);
    nLocal = ncoord*cos(theta)-ecoord*sin(theta);
    eLocal = ncoord*sin(theta)+ecoord*cos(theta);
    fill(eLocal,nLocal,'yellow','FaceAlpha',occMap(iPlot(i))-0.5)
end