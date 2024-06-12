function [localMap,fullMap,measMap,occupiedCells] = moveLocalMap(ego,localMap,fullMap,measMap,occupiedCells,model)
    %Go through the whole local map and move all the elements to a new
    %reference point. Transcribe the local map to the offline full map.
    mapSize = size(localMap.N,2);
    egoN = round(ego.N); 
    egoE = round(ego.E);
    mapCenter = round(mapSize/2);
    refN = localMap.N(mapCenter,mapCenter);
    refE = localMap.E(mapCenter,mapCenter);
    diffN = round((egoN-refN)/model.delta);
    diffE = round((egoE-refE)/model.delta);
    newLocalMap = localMap;
    newOccupiedCells = zeros(mapSize);
    if(abs(diffN)>0.4*model.maxRange || abs(diffE)>0.4*model.maxRange)
        newLocalMap = localMap;
        for i = 1:mapSize
            for j = 1:mapSize
                inMap = i+diffN>1 && i+diffN<mapSize;
                inMap = inMap && j+diffE>1 && j+diffE<mapSize;
                if(inMap)
                    newLocalMap.N(i,j) = localMap.N(i+diffN,j+diffE);
                    newLocalMap.E(i,j) = localMap.E(i+diffN,j+diffE);
                    newLocalMap.occProb(i,j) = localMap.occProb(i+diffN,j+diffE);
                    newLocalMap.Pd(i,j) = localMap.Pd(i+diffN,j+diffE);
                    newLocalMap.alpha(i,j) = localMap.alpha(i+diffN,j+diffE);
                    newLocalMap.beta(i,j) = localMap.beta(i+diffN,j+diffE);
                    newLocalMap.corrFactor(i,j) = localMap.corrFactor(i+diffN,j+diffE);
                    fullMap(150+i+round(refN/model.delta),150+j+round(refE/model.delta))=localMap.occProb(i,j);
                    measMap(150+i+round(refN/model.delta),150+j+round(refE/model.delta))=localMap.alpha(i,j)/localMap.beta(i,j);
                    newOccupiedCells(i,j) = occupiedCells(i+diffN,j+diffE);
                else
                    newLocalMap.alpha(i,j) = 800*model.delta*model.eta;
                    newLocalMap.beta(i,j) = 10*model.eta;
                    newLocalMap.occProb(i,j) = model.initProb;
                    newLocalMap.corrFactor(i,j) = 1;
                    newLocalMap.Pd(i,j) = model.Pd;
                    newLocalMap.N(i,j) = round(i*model.delta+egoN)-round(mapSize/2)*model.delta;
                    newLocalMap.E(i,j) = round(j*model.delta+egoE)-round(mapSize/2)*model.delta;
                    %newLocalMap{i,j}{1} = mapCell;
                end
            end
        end
        localMap = newLocalMap;
        occupiedCells = newOccupiedCells;
    end
end



