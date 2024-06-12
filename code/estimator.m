function [estimates,trajectoryEstimates] = estimator(MBM,model)

trajectoryEstimates = [];

estimates.state = zeros(model.Nx,0);
estimates.cov = zeros(model.Nx,model.Nx,0);

[~,idx] = max(MBM.w);
table_entry = MBM.table(idx,:);
for i = 1:length(table_entry)
    if table_entry(i) > 0 && MBM.track{i}(table_entry(i)).Bern.r >= model.exist_r
        [~,ideath] = max(MBM.track{i}(table_entry(i)).Bern.w_death);
        if ideath == length(MBM.track{i}(table_entry(i)).Bern.w_death)            
            estimates.state = [estimates.state MBM.track{i}(table_entry(i)).Bern.GPState(end).state];
            cov = MBM.track{i}(table_entry(i)).Bern.GPState(end).cov;
            estimates.cov = cat(3,estimates.cov,cov);
        end
        t_death = MBM.track{i}(table_entry(i)).Bern.t_death(ideath);
        tlen = t_death - MBM.track{i}(table_entry(i)).Bern.t_birth + 1;
        

        trajectoryEstimates(end+1,1).t_birth = [MBM.track{i}(table_entry(i)).assocHistory(1).t];
        trajectoryEstimates(end,1).t_death = t_death;
        for k = 1:tlen
            trajectoryEstimates(end,1).xKin{k} = MBM.track{i}(table_entry(i)).Bern.GPState(k).state(1:model.Nkin);
            [extendNE,lowerextendNE,upperextendNE] = getExtendCoordinates(MBM.track{i}(table_entry(i)).Bern.GPState(k),model);
            trajectoryEstimates(end,1).extendNE{k} = extendNE;
            trajectoryEstimates(end,1).lowerextendNE{k} = lowerextendNE;
            trajectoryEstimates(end,1).upperextendNE{k} = upperextendNE;
        end
    end
end  
    
    
end

