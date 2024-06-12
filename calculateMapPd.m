function [localMap,cellsToUpdate] = calculateMapPd(ego,localMap,occupiedCell,model)
    occlInd = find(occupiedCell);
    occlCells = ~isempty(occlInd);
    localMapSize = size(localMap.N,2);
    N = ego.N;
    E = ego.E;
    if(occlCells)
        nOccl = size(occlInd,1);
        radii = zeros(1,nOccl);
        bearingMax = zeros(1,nOccl);
        bearingMin = zeros(1,nOccl);
        for i = 1:nOccl
            %indCell = localMap{occlInd(i)}{1};
            radii(i) = sqrt((N-localMap.N(occlInd(i)))^2+(E-localMap.E(occlInd(i)))^2);
            deltaYmax = localMap.E(occlInd(i))-E+0.5*model.delta;
            deltaXmax = localMap.N(occlInd(i))-N+0.5*model.delta;
            deltaYmin = localMap.E(occlInd(i))-E-0.5*model.delta;
            deltaXmin = localMap.N(occlInd(i))-N-0.5*model.delta;
            [bearingMax(i),bearingMin(i)] = getCellBearings(deltaYmax,deltaYmin,deltaXmax,deltaXmin);
        end

    end
    visibleCells = false(localMapSize);
    corrFactor = localMap.corrFactor;
    Pd = localMap.Pd;
    occProb = localMap.occProb(occlInd);
    r = sqrt((N-localMap.N).^2+(E-localMap.E).^2);
    theta = atan2(localMap.E-E,localMap.N-N);
    corrFactor = 1./r;
    Pd(r>model.maxRange+6) = 0;
    if(occlCells)
        nSectors = 20;
        sectors=linspace(-pi,pi,nSectors);
        for i = 2:nSectors
            maxCond = bearingMax<sectors(i) & bearingMax>sectors(i-1);
            minCond = bearingMin<sectors(i) & bearingMin>sectors(i-1);
            if(i == nSectors)
                maxCond = maxCond | bearingMax>sectors(i);
%             elseif(i == 2) %Handle wrap around in elegant manner, most important issue
%                 temp = bearingMax-2*pi;
%                 maxCond = maxCond(temp>-pi);
            end
            currOccl = find(maxCond | minCond);
            currOccProb = occProb(currOccl);
            currBmax = bearingMax(currOccl);
            currBmin = bearingMin(currOccl);
            currRadii = radii(currOccl);
            currCells = find(theta<sectors(i) & theta>sectors(i-1));
            [~,iSorted] = sort(r(currCells),'ascend');
            nCells = size(currCells,1);
            if(isempty(currOccl))
                Pd(currCells) = model.Pd;
            else
                for j = 1:nCells
                    cellInd = currCells(iSorted(j));
                    rCond = r(cellInd)>currRadii;
                    completelyOccluded = false;
                    if(j>10 && all(rCond) && r(cellInd)<model.maxRange+6)
                       allInd = currCells(iSorted(j-10:j));
                       completelyOccluded = all(Pd(allInd)<0.05);
                    end
                    if(all(rCond) && r(cellInd)<model.maxRange+6 && completelyOccluded)
                        rCondMargin = r(cellInd)>currRadii+10;
                        if(all(rCondMargin))
                            Pd(currCells(iSorted(j:end))) = 0.0001;
                            break
                        end
                    elseif(any(rCond) && r(cellInd)<model.maxRange+6)
                        bearingCond = currBmin<theta(cellInd) & theta(cellInd) < currBmax;
        %                 deltaYmax = localMap.E(cellInd)-E+0.5*model.delta;
        %                 deltaXmax = localMap.N(cellInd)-N+0.5*model.delta;
        %                 deltaYmin = localMap.E(cellInd)-E-0.5*model.delta;
        %                 deltaXmin = localMap.N(cellInd)-N-0.5*model.delta;
        %                [bMax,bMin] = getCellBearings(deltaYmax,deltaYmin,deltaXmin,deltaXmax);
        %                 condmin = (bMin >= currBmin) & (currBmax >= bMin);
        %                 condmax = (bMax >= currBmin) & (currBmax >= bMax);
        % %                 bMinCond = bMin>=bearingMin;
        %                 bMaxCond = bMax<=bearingMax;
        %                 bPartMin = bearingMax >= bMax & bMax>=bearingMin;
        %                 bPartMax = bearingMin >= bMin & bMin >= bearingMax;
        %                 occlCond = rCond & condmin & condmax;
        %                 partOcclCond = rCond & (condmin | condmax);
                        %iOccl = find(occlCond);
                        occlCond = rCond & bearingCond;
                        if(any(occlCond))
                            %ind = currOccl(occlCond);
                            occlProb = currOccProb(occlCond);
                            pNonOccl = 1-occlProb;
                            pNonOccl = prod(pNonOccl);
                            Pd(cellInd) = model.Pd*pNonOccl;
        %                 elseif(any(partOcclCond))
        %                    ind = occlInd(currOccl(occlCond));
        %                    pNonOccl = 1-occProb(ind);
        %                    pNonOccl = prod(pNonOccl);
        %                    Pd(cellInd) = model.Pd*pNonOccl;
                        else
                           Pd(cellInd) = model.Pd;
                        end
                    elseif(r(cellInd)<model.maxRange+6)
                        Pd(cellInd) = model.Pd;
                    else
                        Pd(cellInd) = 0;
                    end
                end
            end
        end
    end
    localMap.corrFactor = corrFactor;
    localMap.Pd = Pd;
    cellsToUpdate = find(Pd>0.1);
end
