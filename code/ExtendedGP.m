classdef ExtendedGP
    % Extended Kalman Filter class    
    properties
        cv
        gp
        noiseStd
        Nx
        Nkin
        Next
        virtual %Boolean for using virtual measurements to constrain extent
        useNI
        sensorRes
    end
    
    methods 
        function obj = ExtendedGP(cvModel, model)
            obj.cv = cvModel;
            obj.gp = model.gp;
            obj.noiseStd = model.sigma_r;
            obj.Next = model.gp.Ntest;
            obj.virtual = model.GPvirtual;
            obj.useNI = 1;
            obj.sensorRes = model.sensorRes;
            if cvModel.NDoF == 2
                obj.Nkin = 4;                
                obj.Nx = 4 + model.gp.Ntest;
            elseif cvModel.NDoF == 3
                obj.Nkin = 6;
                obj.Nx = 6 + model.gp.Ntest;
            end
        end
        
        function F = F(obj, dt)
            % jacobian of prediction function
            F = blkdiag(obj.cv.F(dt), obj.gp.F(dt));
        end
        
        function Q = Q(obj, dt)
            % additive discrete noise covariance
            Q = blkdiag(obj.cv.Q(dt), obj.gp.Q(dt));
        end

        function [gated,precomputedCell] = gateMeasurements(obj,xPred,zPos,PPred,gatingparam)
            if(obj.cv.NDoF == 3)
                extendStart = 7;
                phi = xPred(3);
                kinvec = zeros(2,3);
            end
            if(obj.cv.NDoF == 2)
                extendStart = 5;
                phi = atan2(xPred(4),xPred(3));
            end
            Nz = size(zPos,2);
            gated = zeros(1,Nz);
            %lik = zeros(1,Nz);
            precomputedCell = cell(Nz,1);
            for nz=1:Nz
                diffvector = zPos(:,nz)-xPred(1:2);
                unitVector = diffvector/norm(diffvector);
                globalAngle = atan2(diffvector(2), diffvector(1));
                zAngle = ssa(globalAngle-phi);
                %Derivatives from Özkan & Wahlström (2015)
                dUnitVector = ((diffvector*diffvector')./ norm(diffvector)^3 - eye(2)./norm(diffvector));
                dArgument = 1/ norm(diffvector)^2 * [diffvector(2) -diffvector(1)];

                Kzz = obj.gp.covarianceMatrix(zAngle, zAngle);
                Kzt = obj.gp.covarianceMatrix(zAngle, obj.gp.theta_test);
                dKzt = obj.gp.dCovarianceMatrix(zAngle, obj.gp.theta_test); 
                Ktz = Kzt';
                Hzt = Kzt/obj.gp.Ktt;
                dHzt = dKzt/obj.gp.Ktt;               
                Rzt = Kzz - Kzt*(obj.gp.Ktt\Ktz);
                
                HNE = eye(2) + dUnitVector*(Hzt*xPred(extendStart:end,1))...
                            + unitVector*dArgument*(dHzt*xPred(extendStart:end,1))';
                HPsi = -unitVector*(dHzt*xPred(extendStart:end,1));
                HRadii = unitVector*Hzt; 
                Hnz = [HNE HPsi kinvec HRadii];
                Rnz = unitVector*Rzt*unitVector'+obj.noiseStd^2*eye(2);
                Snz = Hnz*PPred*Hnz' + Rnz;
                
                %Gating                       
                zPosPredBody = HRadii*xPred(extendStart:end,1);
                predmeas =xPred(1:2)+zPosPredBody;
                vnz = zPos(:,nz)-predmeas;
                gated(nz)=gatingparam>vnz'/Snz*vnz;
                precomputedCell{nz}.H = Hnz;
                precomputedCell{nz}.R = Rnz;
                precomputedCell{nz}.predmeas = predmeas;
            end
        end

        function [Angles, indices] = getOccupiedAngle(obj,xPred)
            extentAngles = atan2(xPred(2)+xPred(obj.Nkin+1:end).*sin(obj.gp.theta_test'+xPred(3)),xPred(1)+xPred(obj.Nkin+1:end).*cos(obj.gp.theta_test'+xPred(3)));
            [minAngle,minI] = min(extentAngles);
            [maxAngle,maxI] = max(extentAngles);
            if(abs(maxAngle-minAngle)>pi)
                tempminAngle = maxAngle;
                maxAngle = minAngle+2*pi;
                tempminI = maxI;
                maxI = minI;
                minAngle = tempminAngle;
                minI = tempminI;
            end
            %Should move this logic elsewhere
            Angles = [minAngle maxAngle];
            idx = [minI maxI];
            indices = [minI maxI];
%             for i = 1:2
%                 if(idx(i) == 1)
%                     angleDiff = abs([extentAngles(obj.Next)-extentAngles(idx(i)) extentAngles(idx(i))-extentAngles(idx(i)+1)]);
%                 elseif(idx(i) == obj.Next)
%                     angleDiff = abs([extentAngles(obj.Next-1)-extentAngles(idx(i)) extentAngles(idx(i))-extentAngles(1)]);
%                 else
%                     angleDiff = abs(diff(extentAngles(idx(i)-1:idx(i))));
%                 end
%                 if(obj.sensorRes>min(angleDiff))
%                     indices(i) = 0;
%                 end
%             end
        end

        function[angles,varAngles] = getVisibilityAngle(obj,xPred,cov)
            Hinit = zeros(2, obj.Nx);
            Rinit = zeros(2);
            varAngles = zeros(1,2);
            [angles,indices] = obj.getOccupiedAngle(xPred);
            extentstart = 7;
            for i = 1:2 %min 1 max 2
                if(indices(i) ~= 0)
                    actIndex = indices(i);
                    angle = obj.gp.theta_test(actIndex);
                    globalAngle = ssa(angle+xPred(3));

                    Kzt = obj.gp.covarianceMatrix(angle, obj.gp.theta_test);
                    Hzt = Kzt/obj.gp.Ktt;

                    denominator = (xPred(1)+Hzt*xPred(extentstart:end)*cos(globalAngle)).^2+(xPred(2)+Hzt*xPred(extentstart:end)*sin(globalAngle)).^2;
                    dX = -(xPred(2)+Hzt*xPred(extentstart:end)*sin(globalAngle))/denominator;
                    dY = (xPred(1)+Hzt*xPred(extentstart:end)*cos(globalAngle))/denominator;
                    dExtent = Hzt*(xPred(1)*sin(globalAngle)-xPred(2)*cos(globalAngle))/denominator;
                    dPhi = (Hzt*xPred(extentstart:end)*(xPred(1)*cos(globalAngle)+Hzt*xPred(extentstart:end)+xPred(2)*sin(globalAngle)))/denominator;
                                     
                    Hinit(i,1:3) = [dX dY dPhi];
                    Hinit(i,extentstart:end) = dExtent;

                    
                    Rinit(i,i) = obj.sensorRes.^2;
                    varAngles(i) = Hinit(i,:)*cov*Hinit(i,:)'+Rinit(i,i);
                end
            end
        end

        function [zPosPred,H,R,zPosExtremes] = virtualAngularMeasurements(obj, xPred,zPos,indices)
                %Virtual measurements that bound the extent, compared with
                %actual minimum and maximum angles from the associated
                %measurements. 
                %
                %[~,indices] = obj.checkNegativeInformation(xPred,zPos);
                Hinit = zeros(2, obj.Nx);
                Rinit = zeros(2);
                zPosPredInit = zeros(1,2);
                zPosExtremesInit = zeros(1,2);
                [zPosAngles,zPosRadii] = cart2pol(zPos(1,:),zPos(2,:));
                [~,maxI] = max(zPosAngles);
                [~,minI] = min(zPosAngles);
                condmax =  indices(2) == 0;
                condmin =  indices(1) == 0;
                extentstart = 7;
                for i = 1:2 %min 1 max 2
                    if(indices(i) ~= 0)
                        actIndex = indices(i);%getActualIndex(obj,xPred,i);
                        angle = obj.gp.theta_test(actIndex);
                        globalAngle = ssa(angle+xPred(3));


                        Kzt = obj.gp.covarianceMatrix(angle, obj.gp.theta_test);
                        Hzt = Kzt/obj.gp.Ktt;

                        denominator = (xPred(1)+Hzt*xPred(extentstart:end)*cos(globalAngle)).^2+(xPred(2)+Hzt*xPred(extentstart:end)*sin(globalAngle)).^2;
                        dX = -(xPred(2)+Hzt*xPred(extentstart:end)*sin(globalAngle))/denominator;
                        dY = (xPred(1)+Hzt*xPred(extentstart:end)*cos(globalAngle))/denominator;
                        dExtent = Hzt*(xPred(1)*sin(globalAngle)-xPred(2)*cos(globalAngle))/denominator;
                        dPhi = (Hzt*xPred(extentstart:end)*(xPred(1)*cos(globalAngle)+Hzt*xPred(extentstart:end)+xPred(2)*sin(globalAngle)))/denominator;
                                            
                        zPosPredInit(:,i) = atan2(xPred(2)+xPred(6+actIndex)*sin(globalAngle),xPred(1)+xPred(6+actIndex)*cos(globalAngle));
                        Hinit(i,1:3) = [dX dY dPhi];
                        Hinit(i,extentstart:end) = dExtent;
                        
                        Rinit(i,i) = obj.sensorRes.^2;
                    end
                end
                if(condmax && condmin)
                    zPosPred = []; H = []; R= []; zPosExtremes = [];
                    return;
                elseif(condmin)
                    zPosPred=zPosPredInit(:,2); H=Hinit(2,:); R=Rinit(2,2);  zPosExtremes=zPosAngles(maxI);
                    return;
                elseif(condmax)
                    zPosPred=zPosPredInit(:,1); H=Hinit(1,:); R=Rinit(1,1);  zPosExtremes=zPosAngles(minI);
                    return;
                else
                    zPosPred=zPosPredInit; H=Hinit; R=Rinit;  zPosExtremes=[zPosAngles(minI) zPosAngles(maxI)];
                    return;
                end
        end
        
        function [zPosPred,H,R,zPosExtremes] = virtualMeasurements(obj, xPred,zPos,indices)
                %Virtual measurements that bound the extent, compared with
                %actual minimum and maximum angles from the associated
                %measurements. 
                %
                %[~,indices] = obj.checkNegativeInformation(xPred,zPos);
                Hinit = zeros(4, obj.Nx);
                Rinit = zeros(4);
                zPosPredInit = zeros(2,2);
                zPosExtremesInit = zeros(2,2);
                [zPosAngles,zPosRadii] = cart2pol(zPos(1,:),zPos(2,:));
                [~,maxI] = max(zPosAngles);
                [~,minI] = min(zPosAngles);
                [zPosExtremesInit(1,1),zPosExtremesInit(2,1)] = pol2cart(zPosAngles(minI),zPosRadii(minI));
                [zPosExtremesInit(1,2),zPosExtremesInit(2,2)] = pol2cart(zPosAngles(maxI),zPosRadii(maxI));
                condmax =  indices(2) == 0;
                condmin =  indices(1) == 0;
                for i = 1:2 %min 1 max 2
                    if(indices(i) ~= 0)
                        actIndex = indices(i);%getActualIndex(obj,xPred,i);
                        angle = obj.gp.theta_test(actIndex);
                        globalAngle = mod(angle+xPred(3),2*pi);
                        unitVector = [cos(globalAngle); sin(globalAngle)];
                        dUnitVector = [-sin(globalAngle); cos(globalAngle)];
                        Kzz = obj.gp.covarianceMatrix(angle, angle);
                        Kzt = obj.gp.covarianceMatrix(angle, obj.gp.theta_test);
                        Ktz = Kzt';
                        Hzt = Kzt/obj.gp.Ktt;
                        Rzt = Kzz - Kzt*(obj.gp.Ktt\Ktz);
                        HNE = eye(2);
                        HPsi = dUnitVector*(Hzt*xPred(obj.Nkin+1:end));
                        HRadii = unitVector*Hzt;
                                            
                        zPosPredBody = HRadii*xPred(obj.Nkin+1:end,1);
                        zPosPredInit(:,i) =xPred(1:2)+zPosPredBody;
                        Hinit(2*i-1:2*i,1:3) = [HNE HPsi];
                        Hinit(2*i-1:2*i,obj.Nkin+1:end) = HRadii;
                        
                        Rinit(2*i-1:2*i,2*i-1:2*i) = unitVector*Rzt*unitVector'...
                                                    + 3*(obj.noiseStd^2*eye(2));
                    end
                end
                if(condmax && condmin)
                    zPosPred = []; H = []; R= []; zPosExtremes = [];
                    return;
                elseif(condmin)
                    zPosPred=zPosPredInit(:,2); H=Hinit(3:4,:); R=Rinit(3:4,3:4);  zPosExtremes=zPosExtremesInit(:,2);
                    return;
                elseif(condmax)
                    zPosPred=zPosPredInit(:,1); H=Hinit(1:2,:); R=Rinit(1:2,1:2);  zPosExtremes=zPosExtremesInit(:,1);
                    return;
                else
                    zPosPred=zPosPredInit; H=Hinit; R=Rinit;  zPosExtremes=zPosExtremesInit;
                    return;
                end
        end

        function index = getActualIndex(obj,xPred,whichExtreme)
            extentAngles = atan2(xPred(2)+xPred(obj.Nkin+1:end).*sin(obj.gp.theta_test'+xPred(3)),xPred(1)+xPred(obj.Nkin+1:end).*cos(obj.gp.theta_test'+xPred(3)));
            [minAngle,minI] = min(extentAngles);
            [maxAngle,maxI] = max(extentAngles);
            if(abs(maxAngle-minAngle)>pi)
                minAngle = maxAngle;
                maxAngle = minAngle+2*pi;
                minI = maxI;
                maxI = minI;
            end
            if whichExtreme == 1
                index = minI;
            elseif whichExtreme == 2
               index = maxI;
            end
        end

        function [zPosPred, H, R,gated] = predictedMeasurements(obj, xPred, zPos)
            % measurement function and its jacobian
            if(obj.cv.NDoF == 3)
                extendStart = 7;
                phi = xPred(3);
            end
            if(obj.cv.NDoF == 2)
                extendStart = 5;
                phi = atan2(xPred(4),xPred(3));
            end
            Nz = size(zPos,2);
            zPosPred = zeros(2,Nz);
            H = zeros(2*Nz, obj.Nx);
            R = zeros(2*Nz);
            gated = zeros(1,Nz);
            zAngles = zeros(Nz);
            for nz=1:Nz
                diffvector = zPos(:,nz)-xPred(1:2);
                unitVector = diffvector/norm(diffvector);
                globalAngle = atan2(diffvector(2), diffvector(1));
                zAngle = ssa(globalAngle-phi);
                zAngles(nz) = zAngle;
                %Derivatives from Özkan & Wahlström (2015)
                dUnitVector = ((diffvector*diffvector')./ norm(diffvector)^3 - eye(2)./norm(diffvector));
                dArgument = 1/ norm(diffvector)^2 * [diffvector(2) -diffvector(1)];

                Kzz = obj.gp.covarianceMatrix(zAngle, zAngle);
                Kzt = obj.gp.covarianceMatrix(zAngle, obj.gp.theta_test);
                dKzt = obj.gp.dCovarianceMatrix(zAngle, obj.gp.theta_test); 
                Ktz = Kzt';
                Hzt = Kzt/obj.gp.Ktt;
                dHzt = dKzt/obj.gp.Ktt;               
                Rzt = Kzz - Kzt*(obj.gp.Ktt\Ktz);
                
                HNE = eye(2) + dUnitVector*(Hzt*xPred(extendStart:end,1))...
                            + unitVector*dArgument*(dHzt*xPred(extendStart:end,1))';
                HPsi = -unitVector*(dHzt*xPred(extendStart:end,1));
                HRadii = unitVector*Hzt; 
                H(2*nz-1:2*nz,1:3) = [HNE HPsi];
                H(2*nz-1:2*nz,extendStart:end) = HRadii;
                
                R(2*nz-1:2*nz,2*nz-1:2*nz) = unitVector*Rzt*unitVector'...
                                            + obj.noiseStd^2*eye(2);
                                    
                zPosPredBody = HRadii*xPred(extendStart:end,1);
                zPosPred(:,nz) =xPred(1:2)+zPosPredBody;
            end
        end
        
        function [xPred, PPred] = predict(obj, xUpd, PUpd, dt)
            % returns the predicted mean and covariance for a time step Ts
            Fk = obj.F(dt);
            xPred = Fk*xUpd;
            xPred(3) = ssa(xPred(3));
            PPred = Fk*PUpd*Fk' + obj.Q(dt);
        end

        function [xUpd, PUpd, measLik] = update(obj,xPred,PPred,zPos,model)
            % returns the mean and covariance after conditioning on the
            % measurement, and the innovation and the innovation covariance
            if(obj.observability(xPred,zPos,model.Ts)<6 && model.observatilityCriteria) %Check observability based on current measurements (not used)
                maxIter = 1;
            else
                maxIter = model.maxIterGP;
            end
            best = cell(1,6);
            eps=1000;
            bestEps = 10001;
            i = 1;
            %Correct centroid position based on angle and range of
            %measurements
            %if(obj.useNI == 1)
                [xIter,PPred] = obj.correctCentroid(xPred,PPred,zPos,model);
           % else
                %xIter = xPred;
            %end
            if(isnan(xIter(1)))
                xIter = xPred;
            end
            indicesNI = zeros(1,2);
            useNegativeInformation = (obj.useNI == 1  && obj.virtual && all(xIter == xPred) && size(zPos,2)>2);
            if(useNegativeInformation)
                [useNegativeInformation,indicesNI] = obj.checkNegativeInformation(xPred,zPos);
            end
            if(abs(xPred(6))>pi)
                xPred(6) = 0;
            end
            optimumFound = false;
            %Optimization loop for GN
            while(i<=maxIter && eps>0.001)
                if isempty(zPos)
                    xUpd = xPred;
                    PUpd = PPred;
                    measLik = log(0);
                    break
                else
                    [xUpd,Sk,Hk,Wk,Rk,zPosPred] = obj.singleUpdate(xPred,xIter,PPred,zPos,useNegativeInformation,indicesNI);
                    eps = sum(abs(xIter-xUpd));
                    if(eps<bestEps)
                        if(abs(xUpd(6)) > pi/8 && useNegativeInformation)% Negative information causes "spins" do not use negative information if that is the case.
                            useNegativeInformation = false;
                            i = 1;
                            bestEps = 1001;
                            eps = 1000;
                            xIter = xPred;
                            continue;
                        end
                        bestEps = eps;
                        best = {xUpd,Sk,Hk,Wk,zPosPred,Rk,xIter};
                        optimumFound = true;
                    elseif(isempty(xUpd) || isnan(xUpd(1)))
                        xUpd = xPred;
                        PUpd = PPred;
                        measLik = log(0);
                        break
                    end
                    xIter = xUpd;
                    i = i+1;
                    if(i>maxIter && useNegativeInformation && eps>0.1)% Do not use negative information if optimization did not converge.
                         useNegativeInformation = false;
                         i = 1;
                         eps = 1000;
                         bestEps = 1001;
                         xIter = xPred;
                    end
                end
            end
            %Calculate Covariance and likelihood based on found solution
            if(optimumFound)
                xUpd = best{1}; Sk = best{2}; Hk = best{3}; Wk = best{4}; zPosPred=best{5}; Rk=best{6}; xIter = best{7};
                I = eye(size(PPred));
                PUpd = (I-Wk*Hk)*PPred*(I-Wk*Hk)' + Wk*Rk*Wk';
                PUpd = (PUpd+PUpd')/2;
                Sk = (Sk+Sk')/2;
                if(useNegativeInformation)
                    [vPredMeas, ~, ~,vMeas] = obj.virtualAngularMeasurements(xIter, zPos,indicesNI);
                    reZPos = [reshape(zPos,[],1); vMeas'];%reshape([zPos vMeas],[],1);
                    reZPosPred = [reshape(zPosPred,[],1); vPredMeas'];%reshape([zPosPred vPredMeas],[],1);
                else
                    reZPosPred = reshape(zPosPred,[],1);
                    reZPos = reshape(zPos,[],1);
                end
                try 
                    measLik = log_mvnpdf(reZPos,reZPosPred,Sk);
                catch ME
                    if(strcmp(ME.identifier,'stats:mvnpdf:BadMatrixSigma'))
                        measLik = -1e6;
                        xUpd = xPred;
                        PUpd = PPred;
                    end
                end
            end
        end

        function [useNegativeInformation,indices] = checkNegativeInformation(obj,xPred,zPos)
            useNegativeInformation = true;
            [occupiedAngles,indices] = obj.getOccupiedAngle(xPred);
        end

        function Obsrank = observability(obj,x,zPos,dt)
            [~, Hk, ~] = obj.predictedMeasurements(x, zPos);
            Obs = Hk;
            F = obj.F(dt);
            d=size(F,1);
            for i = 1:d-1
                Obs = [Obs;Hk*F^i ];
            end
            Obsrank = rank(Obs);
        end

        function [xUpd, PUpd,epsvec,stateIndex,lik1] = updateDebug(obj, xPred,xStart, PPred,zPospre,maxIter)
            % returns the mean and covariance after conditioning on the
            % measurement, and the innovation and the innovation covariance
            xIter = xStart;
            best = cell(1,5);
            eps=1000;
            bestEps = 10001;
            lik = 2000;
            bestLik = 2001;
            i = 1;
            epsvec = zeros(maxIter,1);
            factor =1.5;
            damping=0.01;
            lik1 = 10000*ones(maxIter,1);
            stateIndex = zeros(maxIter,obj.Nx);
            gated = obj.gateMeasurements(xIter, zPospre,PPred,13);
            zPos = zPospre(:,logical(gated));
            optimumFound = false;
            if(size(zPos,2)>2)
                xIter = obj.correctCentroid(xIter,zPos);
            end
            initlik = obj.costFunction(xPred,PPred,zPos,xIter);
            while(i<=maxIter && eps>0.001)
                if isempty(zPos)
                    xUpd = xPred;
                    PUpd = PPred;
                    zPosPred = [];
                    Sk = [];
                    vk = [];
                    break
                else
%                     alpha=linspace(0.001,1,25);
%                     likval = zeros(25,1);
%                     [xUpd,Sk,Hk,Wk,Rk,Phat,vk] = LMUpdate(obj,xPred,xIter,PPred,zPos,0);
%                     prevlik=(vk'/Rk*vk+(xPred-xUpd)'/PPred*(xPred-xUpd))/2;
%                     [xUpd,Sk,Hk,Wk,Rk,Phat,vk] = LMUpdate(obj,xPred,xIter,PPred,zPos,damping);
%                     lik=(vk'/Rk*vk+(xPred-xUpd)'/PPred*(xPred-xUpd))/2;                    
%                     [xUpd,Sk,Hk,Wk,Rk,Phat,vk] = LMUpdate(obj,xPred,xIter,PPred,zPos,damping/factor);
%                     factorlik=(vk'/Rk*vk+(xPred-xUpd)'/PPred*(xPred-xUpd))/2;                    
% 
%                     if(lik<=prevlik)
%                         damping=damping/factor;
%                     elseif(lik>prevlik && factorlik<=prevlik)
%                         damping=damping;
%                     else
%                         while(factorlik>prevlik && damping<1e15)
%                             damping=damping*factor;
%                             [xUpd,Sk,Hk,Wk,Rk,Phat,vk] = LMUpdate(obj,xPred,xIter,PPred,zPos,damping);
%                             factorlik=(vk'/Rk*vk+(xPred-xUpd)'/PPred*(xPred-xUpd))/2;
%                         end
%                     end
%                     [xUpd,Sk,Hk,Wk,Rk,Phat,vk] = LMUpdate(obj,xPred,xIter,PPred,zPos,damping);
%                     
%                      deltaX=xUpd-xIter;
%                     for j = 1:25
%                         xAlpha = xIter+alpha(j)*deltaX;
%                         [zPosPred, Hk, Rk] = obj.predictedMeasurements(xAlpha, zPos);
%                         vk = reshape(zPos-zPosPred,[],1);
%                         likval(j) = (vk'/Rk*vk)+(xPred-xAlpha)'/PPred*(xPred-xAlpha);
%                     end
%                     [lik,idx] = min(likval);
%                     xUpd=xIter+alpha(idx)*deltaX;
                    [xUpd,Sk,Hk,Wk,~,zPosPred] = obj.singleUpdate(xPred,xIter,PPred,zPos);

                    lik = obj.costFunction(xPred,PPred,zPos,xIter);
                    if(i>49)
                        'stop'
                    end
                    eps = sum(abs(xIter-xUpd));
                    [~,index] = max(abs(xIter-xUpd));
                    epsvec(i) = eps;
                    stateIndex(i,:)=xUpd;
                    if(eps<bestEps)
                        bestEps = eps;
                        best = {xUpd,Sk,Hk,Wk,zPosPred};
                    optimumFound = true;
                    elseif(isempty(xUpd) || isnan(xUpd(1)))
                        xUpd = xPred;
                        PUpd = PPred;
                        break
                    end
                    lik1(i) = lik;
                    xIter = xUpd;
                    %Iterating between two states without any improvement,
                    %just stop the iteration (Usually means
                    %divergence)
                    if(i>10)
                        diffs = [diff(epsvec); 0];
                        diffs2 = diff(epsvec(1:2:end));
                        if(any(abs(diffs(1:i))<0.00001))% || any(abs(diffs2(1:floor(i/2)))<0.00001))
                            break
                        end
                    end
                    i = i+1;
                end
            end
            if(optimumFound)
                %PPred = Phat;
                 xUpd = best{1}; Sk = best{2}; Hk = best{3}; Wk = best{4}; zPosPred=best{5};
                I = eye(size(PPred));
                PUpd = (I-Wk*Hk)*PPred;%*(I-Wk*Hk)' + Wk*Rk*Wk';
                %PUpd = (I-Wk*Hk)*PPred;
                PUpd = (PUpd+PUpd')/2;
                k = 2*size(zPos,2);
                epsvec(i-1)=1337;
                %lik1 = lik1-(vk'/Sk*vk+log(det(Sk))+k*log(2*pi))/2;
                %lik1 = (vk'/Rk*vk+(xPred-xUpd)'/PPred*(xPred-xUpd))/2;
                %+log(det(Sk))+k*log(2*pi))/2;
                %lik1=observability(obj.F(0.2),Hk);              
            end
        end

        function [xUpd,Sk,Hk,Wk,Rk,zPosPred] = singleUpdate(obj,xPred,xIter,PPred,zPos,useNegativeInformation,indicesNI)
            [zPosPred, Hk, Rk] = obj.predictedMeasurements(xIter, zPos);
            if(useNegativeInformation)
                [vPredMeas, Hvirt, RVirt,vMeas] = obj.virtualAngularMeasurements(xIter, zPos,indicesNI);
                vk = [reshape(zPos-zPosPred,[],1); vMeas'-vPredMeas']; %reshape([zPos vMeas]-[zPosPred vPredMeas],[],1);
                Hk = [Hk; Hvirt];
                Rk = blkdiag(Rk,RVirt);
            else
                vk = reshape(zPos-zPosPred,[],1);
            end
            Sk = Hk*PPred*Hk' + Rk;
            Wk = PPred*Hk'/Sk;
            xUpd = xPred + Wk*(vk-Hk*(xPred-xIter));
        end

        function [xUpd,Sk,Hk,Wk,Rk,Phat,vk] = LMUpdate(obj,xPred,xIter,PPred,zPos,damping)
            [zPosPred, Hk, Rk] = obj.predictedMeasurements(xIter, zPos);
            vk = reshape(zPos-zPosPred,[],1);
            J2k = Hk'/Rk*Hk+inv(PPred);
            Bk = diag(diag(J2k));
            I = eye(size(PPred));

            %Phat = (I-PPred/(PPred+(1/damping)*inv(Bk)))*PPred;
            Phat=inv(inv(PPred)+damping*Bk);
            Sk = Hk*Phat*Hk' + Rk;
            Wk = Phat*Hk'/Sk;
            xUpd = xPred + Wk*(vk-Hk*(xPred-xIter))-damping*(I-Wk*Hk)*Phat*Bk*(xPred-xIter);
            lik=(vk'/Rk*vk+(xPred-xUpd)'/PPred*(xPred-xUpd))/2;
        end

       
        function lik = calculateLogLikelihood(obj,x,P,zPos,lik_matrix)
            nMeas = length(lik_matrix);
            H = zeros(nMeas*2,obj.Nx);
            R = zeros(nMeas*2);
            predMeas = zeros(2,nMeas);
            for i = 1:nMeas
                H(2*i-1:2*i,:) = lik_matrix{i}.H;
                predMeas(:,i) = lik_matrix{i}.predmeas;
                R(2*i-1:2*i,2*i-1:2*i) = lik_matrix{i}.R;
            end
            reZPosPred = reshape(predMeas,[],1);
            reZPos = reshape(zPos,[],1);
            S = H*P*H'+R;
            S = (S+S')/2;
            try 
                lik = log_mvnpdf(reZPos,reZPosPred,S);
            catch ME
                if(strcmp(ME.identifier,'stats:mvnpdf:BadMatrixSigma'))
                    lik = -1e6;
                end
            end
        end
        
        function lik = calculateLogLikelihoodNegative(obj,x,P,zPos,lik_matrix,useNI)
            nMeas = length(lik_matrix);
            H = zeros(nMeas*2,obj.Nx);
            R = zeros(nMeas*2);
            zPosPred = zeros(2,nMeas);
            for i = 1:nMeas
                H(2*i-1:2*i,:) = lik_matrix{i}.H;
                zPosPred(:,i) = lik_matrix{i}.predmeas;
                R(2*i-1:2*i,2*i-1:2*i) = lik_matrix{i}.R;
            end
            useNegativeInformation = (useNI == 1  && obj.virtual && size(zPos,2)>2);
            if(useNegativeInformation)
                [useNegativeInformation,indicesNI] = obj.checkNegativeInformation(x,zPos);
            end
            if(useNegativeInformation)
                [vPredMeas, Hvirt, RVirt,vMeas] = obj.virtualAngularMeasurements(x, zPos, indicesNI);
                reZPos = [reshape(zPos,[],1); vMeas'];%reshape([zPos vMeas],[],1);
                reZPosPred = [reshape(zPosPred,[],1); vPredMeas'];%reshape([zPosPred vPredMeas],[],1);
                Hk = [H; Hvirt];
                Rk = blkdiag(R,RVirt);
                Sk = Hk*P*Hk'+Rk;
            else
                Sk = H*P*H'+R;
                reZPosPred = reshape(zPosPred,[],1);
                reZPos = reshape(zPos,[],1);                 
            end
            Sk = (Sk+Sk')/2;
            try 
                lik = log_mvnpdf(reZPos,reZPosPred,Sk);
            catch ME
                if(strcmp(ME.identifier,'stats:mvnpdf:BadMatrixSigma'))
                    lik = -1e6;
                end
            end
        end

        function lik = measLogLikelihood(obj,x,P,zPos)
            %Calculate likelihood of measurements
            [zPosPred, Hk, Rk] = obj.predictedMeasurements(x, zPos);
            Sk = Hk*P*Hk'+Rk;
            reZPosPred = reshape(zPosPred,[],1);
            reZPos = reshape(zPos,[],1);
            Sk = (Sk+Sk')/2;
            try 
                lik = log_mvnpdf(reZPos,reZPosPred,Sk);
            catch ME
                if(strcmp(ME.identifier,'stats:mvnpdf:BadMatrixSigma'))
                    lik = -1e6;
                end
            end
        end

        function lik = measLogLikelihoodNegative(obj,x,P,zPos,useNI)
            %New method to handle virtual measurements, not done for
            %update, in general needs to be better, like calling it with
            %the whole struct.
            [zPosPred, Hk, Rk] = obj.predictedMeasurements(x, zPos);
            useNegativeInformation = (useNI == 1  && obj.virtual && size(zPos,2)>2);
            if(useNegativeInformation)
                [useNegativeInformation,indicesNI] = obj.checkNegativeInformation(x,zPos);
            end
            if(useNegativeInformation)
                [vPredMeas, Hvirt, RVirt,vMeas] = obj.virtualAngularMeasurements(x, zPos, indicesNI);
                reZPos = [reshape(zPos,[],1); vMeas'];%reshape([zPos vMeas],[],1);
                reZPosPred = [reshape(zPosPred,[],1); vPredMeas'];%reshape([zPosPred vPredMeas],[],1);
                Hk = [Hk; Hvirt];
                Rk = blkdiag(Rk,RVirt);
                Sk = Hk*P*Hk'+Rk;
            else
                Sk = Hk*P*Hk'+Rk;
                reZPosPred = reshape(zPosPred,[],1);
                reZPos = reshape(zPos,[],1);
                    
            end
            Sk = (Sk+Sk')/2;
            try 
                lik = log_mvnpdf(reZPos,reZPosPred,Sk);
            catch ME
                if(strcmp(ME.identifier,'stats:mvnpdf:BadMatrixSigma'))
                    lik = -1e6;
                end
            end
        end

        function lik = costFunction(obj,xPred,PPred,zPos,x)
            [zPosPred, Hk, Rk] = obj.predictedMeasurements(x, zPos);
            useNegativeInformation = (obj.useNI == 1 && size(zPos,2)>2  && obj.virtual);
            if(useNegativeInformation)
                [useNegativeInformation,indicesNI] = obj.checkNegativeInformation(xPred,zPos);
            end
            if(useNegativeInformation)
                [vPredMeas, Hvirt, RVirt,vMeas] = obj.virtualAngularMeasurements(x, zPos,indicesNI);
                reZPos = [reshape(zPos,[],1); vMeas];%reshape([zPos vMeas],[],1);
                reZPosPred = [reshape(zPosPred,[],1); vPredMeas];%reshape([zPosPred vPredMeas],[],1);
                Rk = blkdiag(Rk,RVirt);
            else
                reZPosPred = reshape(zPosPred,[],1);
                reZPos = reshape(zPos,[],1);
            end
            xi = [reshape(reZPosPred,[],1); x];
            preds = [reshape(reZPos,[],1); xPred];
            Qk = blkdiag(Rk,PPred);
             try 
                lik = log_mvnpdf(xi,preds,Qk);
            catch ME
                if(strcmp(ME.identifier,'stats:mvnpdf:BadMatrixSigma'))
                    lik = -1e6;
                end
            end

        end
        function [xIter,PUpd] = correctCentroid(obj,xIter,PPred,zPos,model)
            [centAngle,centDistance]=cart2pol(xIter(1),xIter(2));
            [zAngle,zRadii] = cart2pol(zPos(1,:),zPos(2,:));          
            anglecond = centAngle<min(zAngle) || max(zAngle)<centAngle;
            distCond = centDistance<min(zRadii);
            PUpd = PPred;
            headingCond = angleDistance(xIter(3),atan2(xIter(5),xIter(4)))>pi/3;
            headingCond = headingCond && sqrt(xIter(4)^2+xIter(5)^2)>1;
            angVelCond = false;%abs(xIter(6))>pi/8;
            if(anglecond && size(zPos,2)>2 && model.positionCorrection)
                %radii = centDistance;
                radii = mean(zRadii)+min(xIter(obj.Nkin+1:end));
                if(centAngle > max(zAngle))
                    angle = max(zAngle);
                else
                    angle = min(zAngle);
                end
                angle = mean(zAngle); 
                [x,y] = pol2cart(angle,radii);
                xIter(1) = x;
                xIter(2) = y;
            end
            if(distCond && size(zPos,2)>2 && model.positionCorrection)
                radii = mean(zRadii)+min(xIter(obj.Nkin+1:end));
                %radii = min(zRadii)+min(xIter(obj.Nkin+1:end));
                %angle = atan2(xIter(2),xIter(1));
                angle = mean(zAngle);
                [x,y] = pol2cart(angle,radii);
                xIter(1) = x;
                xIter(2) = y;              
            end
            if((headingCond || angVelCond) && size(zPos,2)>2 && model.headingCorrection)
                xIter(3) = ssa(atan2(xIter(5),xIter(4)));
                %PUpd(:,3) = zeros(1,obj.Nx);
                %PUpd(3,:) = zeros(1,obj.Nx);
                PUpd(3,3) = PUpd(3,3)*10;
                xIter(6) = 0;
                %PUpd(:,6) = zeros(1,obj.Nx);
                %PUpd(6,:) = zeros(1,obj.Nx);
                %PUpd(6,6) = pi/2;
            end
%             if(angVelCond && size(zPos,2)>2 && model.headingCorrection)
%                  %[xIter,PUpd] = obj.checkHeading(xIter,PPred,zPos);
%             end
        end

        function [xUpd,PUpd] = checkHeading(obj,xPred,PPred,zPos)
            %Maximum likelihood check of heading by checking each test
            %angle based on measurement likelihood only. (Not used)
            inds = linspace(1,obj.Next,obj.Next);
            bestI = 1;
            bestLik = -11111;
            xUpd = xPred;
            PUpd = PPred;
            x=xPred;
            P = PPred;
            for i = 1:obj.Next
                ind = 1+mod(inds-2+i,obj.Next);
                x(3) = ssa(xPred(3)+obj.gp.theta_test(i));
%                 x(obj.Nkin+1:end) = xPred(obj.Nkin+ind);
%                 P(obj.Nkin+1:end,:) = PPred(obj.Nkin+ind,:);
%                 P(:,obj.Nkin+1:end) = PPred(:,obj.Nkin+ind);
                lik = obj.measLogLikelihood(x,P,zPos);
                if(lik>bestLik)
                    bestLik = lik;
                    bestI = i;
                end
            end
            ind = 1+mod(inds-2+bestI,obj.Next);
            xUpd(3) = mod(xPred(3)+obj.gp.theta_test(bestI),2*pi);
            PUpd(:,3) = zeros(1,obj.Nx);
            PUpd(3,:) = zeros(1,obj.Nx);
            PUpd(3,3) = PPred(3,3);
%             xUpd(obj.Nkin+1:end) = xPred(obj.Nkin+ind);
%             PUpd(obj.Nkin+1:end,:) = PPred(obj.Nkin+ind,:);
%             PUpd(:,obj.Nkin+1:end) = PPred(:,obj.Nkin+ind);
%             PUpd(obj.Nkin+1:end,obj.Nkin+1:end) =  PPred(obj.Nkin+ind,obj.Nkin+ind);
        end


    end
end