function [PPP,MBM,localMap,occupiedCells,times] = updatePMBM(PPP,MBM,zPos,time,model,filter,localMap,occupiedCells,egoState)
tic
[localMap,cellsToUpdate] = calculateMapPd(egoState,localMap,occupiedCells,model);
times = zeros(7,1);
[PPP,MBM] = occlusionStaticDynamic(MBM,PPP,egoState,localMap,time,filter,model);
%Update map occlusion
if model.occlusion
    [PPP,MBM] = occlusionCorrection(MBM,PPP,time,filter,model,egoState);
end
[localMap] = occlusionDynamicStatic(MBM,time,filter,model,localMap,egoState,cellsToUpdate);
times(1) = toc;
tic
nMeas = size(zPos,2);                      %number of measurements received
used_meas_u = false(nMeas,1);           %measurement indices inside the gate of undetected objects
nu = length(PPP.w);                 %number of mixture components in PPP intensity
gating_matrix_u = false(nMeas,nu);      %gating matrix for PPP components
lik_matrix_u = cell(nMeas,nu);
for i = 1:nu
    %Perform gating for each mixture component in the PPP intensity
    [gating_matrix_u(:,i),lik_matrix_u(:,i)] = filter.gateMeasurements(PPP.GPState(i).state,zPos,PPP.GPState(i).cov,model.gamma);
    used_meas_u = used_meas_u | gating_matrix_u(:,i);
end
n_tt = length(MBM.track);           %number of pre-existing tracks
gating_matrix_d = cell(n_tt,1);     %gating matrix for each track
used_meas_d = false(nMeas,1);           %measurement indices inside the gate of detected objects
lik_matrix_d = cell(n_tt,1);
for i = 1:n_tt
    %number of hypotheses in track i
    num_hypo = length(MBM.track{i});
    %construct gating matrix
    gating_matrix_d{i} = false(nMeas,num_hypo);
    lik_matrix_d{i} = cell(nMeas,num_hypo);
    for j = 1:num_hypo
        currBern = MBM.track{i}(j).Bern;
        %Perform gating for each single object hypothesis
        if currBern.t_death(end) == time && currBern.w_death(end) > 0
            [gating_matrix_d{i}(:,j),lik_matrix_d{i}(:,j)] = filter.gateMeasurements(currBern.GPState(end).state,zPos,currBern.GPState(end).cov,model.gamma);
        end
    end
    used_meas_d = used_meas_d | sum(gating_matrix_d{i},2) >= 1;
end
%Find the cell associated to each measurement
[detectedCell,measuredCells,gating_matrix_m] = findAssocCell(zPos,localMap,model);
times(2) = toc;
tic;
n_cells = size(measuredCells,1);
%measurement indices inside the gate
used_meas = used_meas_d | used_meas_u;
%find indices of measurements inside the gate of undetected
%objects but not detected objects
used_meas_u_not_d = used_meas > used_meas_d;
%Find measurements that are only used by the map
used_meas_m_only = ~used_meas;

sharedMapDetected = any(used_meas_d & gating_matrix_m, 1);
iCellsDetectedTargets = find(sharedMapDetected);

sharedMapUndetected = any(used_meas_u_not_d & gating_matrix_m, 1);
iCellsUndetectedTargets = find(sharedMapUndetected);
iCellsTargets = find((sharedMapDetected | sharedMapUndetected));
iCellsMapOnly = find(~(sharedMapDetected | sharedMapUndetected));

%used measurements by detected targets
W1 = zPos(:,used_meas);
gating_matrix_u1 = gating_matrix_u(used_meas,:);
gating_matrix_d = cellfun(@(x) x(used_meas,:), gating_matrix_d, 'UniformOutput',false);
lik_matrix_u1 = lik_matrix_u(used_meas_d,:);
lik_matrix_d = cellfun(@(x) x(used_meas,:), lik_matrix_d, 'UniformOutput',false);

%Data association using stochastic optimisation
J = length(MBM.w);
if(J>0)
    wAssoc = [];
    Nj = zeros(J,1);
    P = cell(J,1);
    for j = 1:J
        if isempty(W1)
            %Misdetection
            track_indices = find(MBM.table(j,:)>0);
            nj = length(track_indices);
            P{j}{1} = cell(nj,1);
            lik = 0;
            for i = 1:nj
                %Create new Bernoulli component due to missed detection
                track_miss = MBM.track{track_indices(i)}(MBM.table(j,track_indices(i)));
                %Compute the misdetection likelihood
                [~,lik_miss] = misdetectionBern(track_miss.Bern,model);
                lik = lik + lik_miss;
            end
        else
            %Find the most likely global hypotheses using Stochastic
            %Optimization
            [P{j},lik] = ObjectsSO(MBM,j,PPP,W1,gating_matrix_d,lik_matrix_d,gating_matrix_u(used_meas,:),lik_matrix_u(used_meas,:),measuredCells(iCellsTargets,:)...
                ,gating_matrix_m(used_meas,iCellsTargets),localMap,model,filter);
        end
        Nj(j) = length(lik);
        wAssoc = [wAssoc;lik+MBM.w(j)];
    end
else
   [P{1},lik] = ObjectsSO(MBM,0,PPP,W1,gating_matrix_d,lik_matrix_d,gating_matrix_u(used_meas,:),lik_matrix_u(used_meas,:),measuredCells(iCellsTargets,:)...
                ,gating_matrix_m(used_meas,iCellsTargets),localMap,model,filter);
   wAssoc = lik;
   Nj = length(lik);
end
%Normalise and prune components with low weights
[wAssoc,~] = normalizeLogWeights(wAssoc);
threshold = log(model.threshold_w);
idx_keep = wAssoc >= threshold;
wAssoc = wAssoc(idx_keep);
if length(wAssoc) == 1; wAssoc = 0; end
times(3) = toc;
tic;
%Remove measurement partitions that correspond to low weight global hypotheses
true_meas_indices = find(used_meas==1);
nMeas = size(zPos,2);
idx_0 = 0;
nMBs = length(P);
for j = 1:nMBs
    idx_j = idx_keep(idx_0+1:idx_0+Nj(j));
    P{j} = P{j}(idx_j);
    %Convert set representation to boolean vector representation
    if ~isempty(P{j})
        for i = 1:length(P{j})
            for k = 1:length(P{j}{i})
                P{j}{i}{k} = ismember(1:nMeas,true_meas_indices(P{j}{i}{k}));
            end
        end
    end
    idx_0 = idx_0 + Nj(j);
end

%For each single target hypothesis find its corresponding measurement cells
meas_track = cell(n_tt,1);
for i = 1:n_tt
    meas_track{i} = cell(length(MBM.track{i}),1);
end
meas_mapCell = cell(n_cells,1);

n_mbCells = size(iCellsTargets,2);
meas_newtracks = false(1,nMeas);
idx_new = 0;
for j = 1:nMBs
    if ~isempty(P{j})
        if(J>0)
            track_indices = find(MBM.table(j,:)>0);
            nj = length(track_indices);
        else
            nj = 0;
        end
        n_Bern = nj+n_mbCells;
        for i = 1:length(P{j})
            for k = 1:nj
                if isempty(meas_track{track_indices(k)}{MBM.table(j,track_indices(k))})
                    meas_track{track_indices(k)}{MBM.table(j,track_indices(k))} = ...
                        [meas_track{track_indices(k)}{MBM.table(j,track_indices(k))};P{j}{i}{k}];
                elseif ~ismember(P{j}{i}{k},meas_track{track_indices(k)}{MBM.table(j,track_indices(k))},'rows')
                    meas_track{track_indices(k)}{MBM.table(j,track_indices(k))} = ...
                        [meas_track{track_indices(k)}{MBM.table(j,track_indices(k))};P{j}{i}{k}];
                end
            end
            for k = 1+nj:n_Bern
                iCell = iCellsTargets(k-nj);
                if isempty(meas_mapCell{iCell})
                    meas_mapCell{iCell} = [meas_mapCell{iCell};P{j}{i}{k}];
                elseif ~ismember(P{j}{i}{k},meas_mapCell{iCell},'rows')
                    meas_mapCell{iCell} = [meas_mapCell{iCell};P{j}{i}{k}];
                end
            end
            num_mea_cell = length(P{j}{i});
            if num_mea_cell > n_Bern
                for k = 1:num_mea_cell-n_Bern
                    if ~ismember(P{j}{i}{k+n_Bern},meas_newtracks,'rows')
                        idx_new = idx_new + 1;
                        meas_track{n_tt+idx_new,1}{1} = P{j}{i}{k+n_Bern};
                        meas_newtracks = [meas_newtracks;P{j}{i}{k+n_Bern}];
                    end
                end
            end
        end
    end
end

%Make sure boolean representation
for i = 1:n_tt
    for j = 1:length(meas_track{i})
        meas_track{i}{j} = logical(meas_track{i}{j});
    end
end
meas_newtracks = logical(meas_newtracks);

%Construct new global hypotheses look-up table
n_tt_upd = length(meas_track);
table = zeros(length(wAssoc),n_tt_upd);
idx = 0;
for j = 1:nMBs
    if ~isempty(P{j})
        if(J>0)
            track_indices = find(MBM.table(j,:)>0);
            nj = length(track_indices);
        else
            nj = 0;
        end
        for i = 1:length(P{j})
            idx = idx + 1;
            for k = 1:nj
                [~,table(idx,track_indices(k))] = ...
                    ismember(P{j}{i}{k},meas_track{track_indices(k)}{MBM.table(j,track_indices(k))},'rows');
                if MBM.table(j,track_indices(k)) > 1
                    for p = 1:MBM.table(j,track_indices(k))-1
                        table(idx,track_indices(k)) = table(idx,track_indices(k)) + size(meas_track{track_indices(k)}{p},1);
                    end
                end
            end
            num_mea_cell = length(P{j}{i});
            n_Bern = nj+n_mbCells;
            if num_mea_cell > n_Bern
                for k = 1:num_mea_cell-n_Bern
                    [~,table_idx] = ismember(P{j}{i}{k+n_Bern},meas_newtracks,'rows');
                    table(idx,table_idx-1+n_tt) = 1;
                end
            end
        end
    end
end
times(4) = toc;
tic;
indices = 1:size(zPos,2);

%Update tracks
tracks = cell(n_tt_upd,1);
for i = 1:n_tt
    idx = 0;
    for j = 1:length(meas_track{i})
        for k = 1:size(meas_track{i}{j},1)
            if (k>1 && size(meas_track{i}{j},2) == 0)
                break;
            else
                idx = idx + 1;
                tracks{i}(idx,1) = MBM.track{i}(j); 
                if any(meas_track{i}{j}(k,:))
                    [Bern,lik] = detectionBern(MBM.track{i}(j).Bern,zPos(:,meas_track{i}{j}(k,:)),model,filter);
                else
                    [Bern,lik] = misdetectionBern(MBM.track{i}(j).Bern,model);
                end
                tracks{i}(idx,1).Bern = Bern;
                tracks{i}(idx,1).lik = tracks{i}(idx,1).lik + lik;
                len = length(tracks{i}(idx,1).assocHistory);
                tracks{i}(idx,1).assocHistory(len+1,1).t = time;
                tracks{i}(idx,1).assocHistory(len+1,1).meas = indices(meas_track{i}{j}(k,:));
            end
        end
    end
end
if n_tt_upd > n_tt
    for i = 1:n_tt_upd - n_tt
        in_gate = sum(gating_matrix_u(meas_track{n_tt+i}{1},:),1)>=1;
        if any(in_gate)
            [Bern,lik] = detectionPPP(zPos(:,meas_track{n_tt+i}{1}),PPP,in_gate,model,filter);
            Bern.t_birth = time;
            Bern.t_death = time;
            Bern.w_death = 1;
        else
            Bern.r = 0; 
            Bern.GPState = struct('state',zeros(1,filter.Nx),'cov',zeros(filter.Nx,filter.Nx),'alpha',0,'beta',1); 
            lik = [];
        end
        tracks{n_tt+i,1} = struct('Bern',Bern,'lik',lik,'assocHistory',[]);
        tracks{n_tt+i,1}.assocHistory(1).t = time;
        tracks{n_tt+i,1}.assocHistory(1).meas = indices(meas_track{n_tt+i}{1});
    end
end
times(5) = toc;
tic;
%Append new tracks
% [tracks,table,wAssoc] = newObjectsSO(tracks,table,wAssoc,PPP,zPos,gating_matrix_u,lik_matrix_u,used_meas_u_not_d,time,measuredCells(iCellsUndetectedTargets,:)...
%             ,gating_matrix_m(:,iCellsUndetectedTargets),localMap,model,filter);
%Need to populate meas_mapCell from three sources and have weights for
%different hypotheses.

times(6) = toc;
for i = 1:size(cellsToUpdate,1)
    ind = cellsToUpdate(i,:);
    if(detectedCell(ind(1),ind(2))>0)
        iCell = find(all(ind == measuredCells,2));
        if(any(iCell==iCellsMapOnly))
            occupiedCells(ind(1),ind(2)) = 1;
            [localMap{ind(1),ind(2)}{1},lik] = updateDetectedCell(localMap{ind(1),ind(2)}{1},detectedCell(ind(1),ind(2)),model); %Simple rules
        else
            occupiedCells(ind(1),ind(2)) = 1;
            [localMap{ind(1),ind(2)}{1},lik] = updateDetectedCellMerge(localMap{ind(1),ind(2)}{1},meas_mapCell{iCell},wAssoc,model); %Simple rules
        end
    else
        [localMap{ind(1),ind(2)}{1},lik] = updateUndetectedCell(localMap{ind(1),ind(2)}{1},model); %Simple rules
    end
end
tic;
%Remove Bernoulli components with low existence probability
n_tt = length(tracks);
for i = 1:n_tt
    %Find all Bernoulli components needed to be pruned
    idx = arrayfun(@(x) x.Bern.r < model.threshold_r, tracks{i});
    %Prune these Bernoulli components
    tracks{i} = tracks{i}(~idx);
    idx = find(idx);
    %Update hypothesis table, if a Bernoulli component is
    %pruned, set its corresponding entry to zero
    for j = 1:length(idx)
        temp = table(:,i);
        temp(temp==idx(j)) = 0;
        table(:,i) = temp;
    end
end

%Remove unused tracks
idx_empty = cellfun('isempty',tracks);
table = table(:,~idx_empty);
tracks = tracks(~idx_empty);

%Remove tracks that contains only null single object hypotheses
idx_keep = sum(table,1) > 0;
table = table(:,idx_keep);
tracks = tracks(idx_keep);
if isempty(table)
    wAssoc = [];
end

%Re-index hypothesis table
n_tt = length(tracks);
for i = 1:n_tt
    idx = table(:,i) > 0;
    [~,~,table(idx,i)] = unique(table(idx,i),'rows','stable');
end

%Merge duplicate hypothesis table rows
if ~isempty(table)
    [ht,~,ic] = unique(table,'rows','stable');
    if(size(ht,1)~=size(table,1))
        %There are duplicate entries
        w = zeros(size(ht,1),1);
        for i = 1:size(ht,1)
            indices_dupli = (ic==i);
            [~,w(i)] = normalizeLogWeights(wAssoc(indices_dupli));
        end
        table = ht;
        wAssoc = w;
    end
end

%%%%%%%%%%%
%Merge similar Bernoulli components in the track (Not used but makes
%runtime faster)

% n_tt = length(tracks);
% for i = 1:n_tt
%     nb = length(tracks{i});
%     if nb > 1
%         %find the highest weight MB that contains this track
%         idx_mb = find(table(:,i) > 0);
%         [~,idx] = max(wAssoc(idx_mb));
%         idx_b = table(idx_mb(idx),i);
%         I = idx_b;
%         for j = 1:nb
%             if j ~= idx_b && tracks{i}(idx_b).Bern.r == tracks{i}(j).Bern.r...
%                     && GP_KLdiv(tracks{i}(idx_b).Bern.GPState(end),tracks{i}(j).Bern.GPState(end)) < model.merge %Note similarity on final is used
%                 I = [I j];
%             end
%         end
%         len = length(I);
%         if len > 1
%             Bern_temp = [tracks{i}(I).Bern];
%             w_temp = zeros(len,1);
%             GPtomerge = [];
%             for l = 1:len
%                 GPtomerge = [GPtomerge; Bern_temp(l).GPState(end)];
%                 [~,w_temp(l)] = normalizeLogWeights(wAssoc(table(:,i)==I(l)));
%             end
%             %Assume the same weight
%             [~,GP_new] = GGP_merge(w_temp,GPtomerge,model);
%             tracks{i}(idx_b).Bern.GPState(end) = GP_new;
%             idx_remove = setdiff(I,idx_b);
%             tracks{i}(idx_remove) = [];
%             table(ismember(table(:,i),I),i) = idx_b;
%         end
%     end
% end
% 
% %Re-index hypothesis table
% n_tt = length(tracks);
% for i = 1:n_tt
%     idx = table(:,i)>0;
%     [~,~,table(idx,i)] = unique(table(idx,i),'rows','stable');
% %     if any(idx)
% %         table(idx,i) = 0;
% %         table(~idx,i) = table(~idx,i) - 1;
% %     end
% end
% 
% %Merge duplicate hypothesis table rows
% if ~isempty(table)
%     [ht,~,ic] = unique(table,'rows');
%     if(size(ht,1)~=size(table,1))
%         %There are duplicate entries
%         w = zeros(size(ht,1),1);
%         for i = 1:size(ht,1)
%             indices_dupli = (ic==i);
%             [~,w(i)] = normalizeLogWeights(wAssoc(indices_dupli));
%         end
%         table = ht;
%         wAssoc = w;
%     end
% end
% 
% if length(wAssoc) == 1; wAssoc = 0; end

%%%%%%%%%%%%

%Assign updated value
MBM.table = table;
MBM.track = tracks;
MBM.w = wAssoc;

%PPP misdetection update
PPP = misdetectionPPP(PPP,model);
%Prune PPP components
idx_keep = PPP.w > log(model.threshold_u);
PPP.w = PPP.w(idx_keep);
PPP.GPState = PPP.GPState(idx_keep);
%times(7) = toc;

end

