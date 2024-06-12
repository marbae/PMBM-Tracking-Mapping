function [] = plotShipLocal(ship,sizeScatter,color,egoState,nt)
    egoHeading = deg2rad(egoState.heading);
    egoN = egoState.N;
    egoE = egoState.E;
    xKin = ship.xKin(:,nt);
    extendNE = ship.extendNE(:,:,nt);
    
    xKin(1) = xKin(1)-egoN;
    xKin(2) = xKin(2)-egoE;
    N = xKin(1)*cos(egoHeading)-xKin(2)*sin(egoHeading);
    E = xKin(1)*sin(egoHeading)+xKin(2)*cos(egoHeading);
    psi = xKin(3)-egoHeading;
    vN = xKin(4)*cos(egoHeading)-xKin(5)*sin(egoHeading);
    vE = xKin(4)*sin(egoHeading)+xKin(5)*cos(egoHeading);
    extendNE(1,:) = extendNE(1,:)-egoN;
    extendNE(2,:) = extendNE(2,:)-egoE;
    extendN = extendNE(1,:)*cos(egoHeading)-extendNE(2,:)*sin(egoHeading);
    extendE = extendNE(1,:)*sin(egoHeading)+extendNE(2,:)*cos(egoHeading);
    plot(extendE,extendN,'color',color)
    scatter(E,N,sizeScatter,color,'filled')
    quiver(E,N,vE,vN,'--','color',color)
    quiver(E,N,sin(psi),cos(psi),'color',color)
end



    