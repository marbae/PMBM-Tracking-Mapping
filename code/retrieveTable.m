function [GOSPAtable,metricsTable] = retrieveTable(GOSPA,metrics,nTargets)
    numMC = size(GOSPA,3);
    GOSPAtable = zeros(4,numMC);
    metricsTable = zeros(7,numMC);
    for i = 1:numMC
        GOSPAtable(:,i) = mean(GOSPA(:,:,i),1);
        GOSPAtable(1,i) = mean(GOSPA((GOSPA(:,1,i)>0),1,i),1);
        GOSPAtable(2,i) = mean(GOSPA((GOSPA(:,2,i)>0),2,i),1);
        for iTar = 1:nTargets
            IOU = [];
            RMSE = [];
            RMSEvel = [];
            headingErr = [];
            NEES = [];
            NEESposvel = [];
            NEESheading = [];
            for j = 1:length(metrics)
                if(~isempty(metrics{j,i}))
                    whichTarget = metrics{j,i}.iGT==iTar;
                    if(any(whichTarget))
                        IOU = [IOU metrics{j,i}.IOU(whichTarget)];
                        RMSE = [RMSE metrics{j,i}.RMSE(whichTarget)];
                        headingErr = [headingErr metrics{j,i}.headingErr(whichTarget)];
                        RMSEvel = [RMSEvel metrics{j,i}.RMSEvel(whichTarget)];
                        NEES = [NEES log(metrics{j,i}.NEES(whichTarget))];
%                         NEESposvel = [NEESposvel log(metrics{j,i}.NEESposvel(whichTarget))];
%                         NEESheading = [NEESheading log(metrics{j,i}.NEESheading(whichTarget))];
                    end
                end
            end
            metricsTable(:,i) = metricsTable(:,i) + [mean(IOU); mean(RMSE); mean(RMSEvel(1:end-1)); mean(headingErr); mean(NEES(1:end-1)); mean(NEESposvel(1:end-1)); mean(NEESheading(1:end-1))];
        end
        metricsTable(:,i) = metricsTable(:,i)./nTargets;
    end
end

