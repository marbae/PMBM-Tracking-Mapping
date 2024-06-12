function lik = GPcostFunction(gp,xPred,PPred,zPos,x)
    noiseStd = 0.1;
    Nz = size(zPos,2);
    zPosPred = zeros(2,Nz);
    Nx = size(xPred,1);
    H = zeros(2*Nz, Nx);
    R = zeros(2*Nz);
    extendStart=7;
    phi = x(3);
    for nz=1:Nz
        diffvector = zPos(:,nz)-x(1:2);
        unitVector = diffvector/norm(diffvector);
        globalAngle = atan2(diffvector(2), diffvector(1));
        zAngle = mod(globalAngle-phi,2*pi);
        %Derivatives from Özkan & Wahlström (2015)
        dUnitVector = ((diffvector*diffvector')./ norm(diffvector)^3 - eye(2)./norm(diffvector));
        dArgument = 1/ norm(diffvector)^2 * [diffvector(2) -diffvector(1)];

        Kzz = gp.covarianceMatrix(zAngle, zAngle);
        Kzt = gp.covarianceMatrix(zAngle, gp.theta_test);
        dKzt = gp.dCovarianceMatrix(zAngle, gp.theta_test); 
        Ktz = Kzt';
        Hzt = Kzt/gp.Ktt;
        dHzt = dKzt/gp.Ktt;               
        Rzt = Kzz - Kzt*(gp.Ktt\Ktz);
        
        HNE = eye(2) + dUnitVector*(Hzt*x(extendStart:end,1))...
                    + unitVector*dArgument*(dHzt*x(extendStart:end,1))';
        HPsi = -unitVector*(dHzt*x(extendStart:end,1));
        HRadii = unitVector*Hzt; 
        H(2*nz-1:2*nz,1:3) = [HNE HPsi];
        H(2*nz-1:2*nz,extendStart:end) = HRadii;
        
        R(2*nz-1:2*nz,2*nz-1:2*nz) = unitVector*Rzt*unitVector'...
                                    + noiseStd^2*eye(2);
                            
        zPosPredBody = HRadii*x(extendStart:end,1);
        zPosPred(:,nz) =x(1:2)+zPosPredBody;
    end
    vk = reshape(zPos-zPosPred,[],1);
    Qk = blkdiag(R,PPred);
    xk = xPred-x;
    v = [vk; xk];
    lik = v'/Qk*v;
end