function [bearingMax,bearingMin] = getCellBearings(yMax,yMin,xMax,xMin)
    if(xMax > 0 && yMax>0)
        bearingMax = atan2(yMax,xMin);
        bearingMin = atan2(yMin,xMax);
    elseif(xMax<0 && yMax>0)
        bearingMax = atan2(yMin,xMin);
        bearingMin = atan2(yMax,xMax);
    elseif(xMax<0 && yMax<0)
        bearingMax = atan2(yMin,xMax);
        bearingMin = atan2(yMax,xMin);
    elseif(xMax>0 && yMax<0)
        bearingMax = atan2(yMax,xMax);
        bearingMin = atan2(yMin,xMin);
    end

    if(abs(bearingMax-bearingMin)>pi && bearingMax<0)
        bearingMax = bearingMax+2*pi;
    end