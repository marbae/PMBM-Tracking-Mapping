function [] = plotShip(xKin,extendNE,sizeScatter,color)
    N = xKin(1);
    E = xKin(2);
    psi = xKin(3);
    vN = xKin(4);
    vE = xKin(5);
    extendN = extendNE(1,:);
    extendE = extendNE(2,:);
    plot(extendE,extendN,'color',color)
    scatter(E,N,sizeScatter,color,'filled')
    quiver(E,N,vE,vN,'--','color',color)
    quiver(E,N,sin(psi),cos(psi),'color',color)
end