function [] = plotEstimatedShip(xKin,extendNE,extendNEUpper,extendNELower,sizeScatter,color)
color_fill = [245,222,179]/255; %wheat    
N = xKin(1);
    E = xKin(2);
    psi = xKin(3);
    vN = xKin(4);
    vE = xKin(5);
    extendN = extendNE(1,:);
    extendE = extendNE(2,:);
    extendNUpper = extendNEUpper(1,:);
    extendEUpper = extendNEUpper(2,:);
    extendNLower = extendNELower(1,:);
    extendELower = extendNELower(2,:);
    fill(extendEUpper,extendNUpper,color_fill);
    fill(extendELower,extendNLower,[1,1,1]);
    plot(extendE,extendN,'color',color)
    scatter(E,N,sizeScatter,color,'filled')
    quiver(E,N,vE,vN,'--','color',color)
    quiver(E,N,sin(psi),cos(psi),'color',color)
end