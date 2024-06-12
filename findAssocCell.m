function [detectedCell,measuredCells,gating_matrix_m] = findAssocCell(zPos,localMap,model)
    %Returns gating matrix and a list of cells to update
    mapSize = size(localMap.N,1);
    detectedCell = zeros(mapSize);
    nMeas = size(zPos,2);
    measToMap = zeros(2,nMeas);
    %refCell = localMap{round(mapSize/2),round(mapSize/2)}{1};
    refPoint = [localMap.N(round(mapSize/2),round(mapSize/2)); localMap.E(round(mapSize/2),round(mapSize/2))];
    refZ = zPos-refPoint;
    refZ = round(refZ/model.delta);
    for i = 1:nMeas
        r = sqrt(refZ(1,i)^2+refZ(2,i)^2);
        %if(r<model.maxRange/model.delta)
            indexN = round(mapSize/2)+refZ(1,i);
            indexE = round(mapSize/2)+refZ(2,i);
            measToMap(:,i) = [indexN indexE];
            detectedCell(indexN,indexE) = detectedCell(indexN,indexE)+1;
            
        %end
    end
    measuredCells = unique(measToMap','rows');
    nCells = size(measuredCells,1);
    gating_matrix_m = false(nMeas,nCells);
    for i = 1:size(measuredCells,1)
        gating_matrix_m(:,i) = all(measuredCells(i,:) == measToMap',2);
    end
end
    