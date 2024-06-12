function [] = plotMeasurements(Z,sizeScatter,color)
    
    if ~isempty(Z) && size(Z,1)>=2
        scatter(Z(2,:),Z(1,:),sizeScatter,color,'LineWidth',1.5)
    elseif ~isempty(Z)
        scatter(Z(:,2),Z(:,1),sizeScatter,color)
    end
end