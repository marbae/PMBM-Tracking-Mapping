function [partitions,lik] = ObjectsSO(MBM,nh,PPP,W,gating_matrix_d,lik_matrix_d,...
    gating_matrix_u,lik_matrix_u,measuredCells,gating_matrix_m,localMap,model,filter)

m = size(W,2);

%Pre-store the results of detectionPPP update to avoid duplicate computation
%%%%%%%%%%%%%

Lstore = -inf;
Mstore = false(1,m);

%%%%%%%%%%%%%

if(nh>0)
    track_indices = find(MBM.table(nh,:)>0);     %extract the jth global hypothesis
    n_tt = length(track_indices);               %number of tracks in the jth global hypothesis
    gating_matrix_dj = zeros(m,n_tt);           %extract gating matrix for each track
else
    n_tt = 0;
end
lik_matrix_dj = cell(m,n_tt);
for i = 1:n_tt
    gating_matrix_dj(:,i) = gating_matrix_d{track_indices(i)}(:,MBM.table(nh,track_indices(i)));
    lik_matrix_dj(:,i) = lik_matrix_d{track_indices(i)}(:,MBM.table(nh,track_indices(i)));
    tracks(i,1) = MBM.track{track_indices(i)}(MBM.table(nh,track_indices(i)));
end

n_cells = size(measuredCells,1);        %Number of map cells associated to measurements

%Pre-computing likelihood of new measurement being a cell of its own
Lc = zeros(m,1);
for i = 1:m
    if any(gating_matrix_u(i,:))
        Lc(i) = likPPP(W(:,i),lik_matrix_u(i,:),PPP,gating_matrix_u(i,:),model,filter);
    else
        Lc(i) = log(model.lambda_fa);
    end
end

%Pre-computing misdetection likelihood
Lmiss = zeros(n_tt+n_cells,1);
for i = 1:n_tt
    [~,Lmiss(i)] = misdetectionBern(tracks(i,1).Bern,model);
end
for i = n_tt+1:n_tt+n_cells
    ind = measuredCells(i-n_tt,:);
    currCell = localMap{ind(1),ind(2)}{1};
    [~,Lmiss(i)] = updateUndetectedCell(currCell,model);
    Lmiss(i) = log(Lmiss(i));
end

%Maximum likelihood initialisation
nu = length(PPP.w);

meas_cell = cell(n_tt+n_cells+nu,1);
ini_likelihood = -inf(n_tt+n_cells+1,m);
for i = 1:n_tt
    log_pD_r = log(model.Pd) + log(tracks(i,1).Bern.r);
    ghat = log(tracks(i,1).Bern.GPState(end).alpha/tracks(i,1).Bern.GPState(end).beta);
    for j = 1:m
        if gating_matrix_dj(j,i)
            ini_likelihood(i,j) = log_pD_r + ghat + filter.measLogLikelihood(tracks(i,1).Bern.GPState(end).state,tracks(i,1).Bern.GPState(end).cov,W(:,j));
        end
    end
end

%Map init
for i = 1:n_cells
    ind = measuredCells(i,:);
    currCell = localMap{ind(1),ind(2)}{1};
    log_pD_r = log(currCell.Pd)+log(currCell.occProb);
    ghat = log((currCell.alpha*currCell.corrFactor)/currCell.beta);
    for j = 1:m
        if gating_matrix_m(j,i)
            ini_likelihood(n_tt+i,j) = log_pD_r+ghat+log(1/model.delta^2);
        end
    end
end

ini_likelihood(n_tt+n_cells+1,:) = Lc';

l = 1;
for j = 1:m
    [~,idx] = max(ini_likelihood(:,j));
    if idx <= n_tt
        meas_cell{idx} = [meas_cell{idx} j];
    elseif idx > n_tt && idx<=n_tt+n_cells
        meas_cell{idx} = [meas_cell{idx} j];
    else
        idx_u = find(gating_matrix_u(j,:)==true,1);
        if ~isempty(idx_u)
            meas_cell{n_tt+n_cells+idx_u,1} = [meas_cell{n_tt+n_cells+idx_u,1} j];
        else
            meas_cell{n_tt+n_cells+nu+l,1} = j;
            l = l + 1;
        end
    end
end

idx = ~cellfun('isempty',meas_cell);
idx(1:n_tt+n_cells) = true;
meas_cell =  meas_cell(idx);

N = length(meas_cell);
Lp = zeros(N,1);
for i = 1:n_tt
    if isempty(meas_cell{i})
        Lp(i) = Lmiss(i);
    else
        Lp(i) = likBern(tracks(i,1).Bern,W(:,meas_cell{i}),lik_matrix_dj(meas_cell{i},i),model,filter);
    end
end
for i = n_tt+1:n_cells+n_tt
    if isempty(meas_cell{i})
        Lp(i) = Lmiss(i);
    else
        ind = measuredCells(i-n_tt,:);
        Lp(i) = likMapCell(localMap{ind(1),ind(2)}{1},W(:,meas_cell{i}),model);
    end
end
%Lp(n_tt+1) = likMapCell(localMap,meas_cell{n_tt+1},model); %Design this in a smart and modular way
for i = n_cells+n_tt+1:N
    in_gate = sum(gating_matrix_u(meas_cell{i},:),1)>=1;
    if any(in_gate)
        Lp(i) = likPPP(W(:,meas_cell{i}),lik_matrix_u(meas_cell{i},:),PPP,in_gate,model,filter);
    else
        Lp(i) = log(model.lambda_fa);
    end
end


T = model.num_iterations*(m+n_tt+n_cells);     %number of iterations
partitions = cell(T,1);
partitions{1} = meas_cell;
lik = zeros(T,1);
lik(1) = sum(Lp);

max_repetition = max(20,ceil(model.max_repetition*(m+n_tt+n_cells)/2));
num_repetition = 0;
%Stochastic optimisation
for t = 1:T
    %randomly select a measurement
    mea_selected = randi(m,1);
    %find the corresponding cell index
    Loc = cellfun(@(x) x==mea_selected,meas_cell,'Un',0);
    cell_idx = find(cellfun(@(x) any(x(:)),Loc));
    N = length(meas_cell);
    
    mea_after_move = meas_cell{cell_idx}(meas_cell{cell_idx}~=mea_selected);
    Wp = -inf(2*N+4,1);
    Lp_original = Lp;
    meas_cell_original = meas_cell;
    
    %selected measurement moved to an existing cell
    Lp_move1 = -inf(N,1);
    if cell_idx <= n_tt+n_cells %selected measurement from a target cell. These if statements should include the map
        if isempty(mea_after_move)
            Lp_move2 = Lmiss(cell_idx);
        else
            if cell_idx <= n_tt
                Lp_move2 = likBern(tracks(cell_idx).Bern,W(:,mea_after_move),lik_matrix_dj(mea_after_move,cell_idx),model,filter);
            else
                ind = measuredCells(cell_idx-n_tt,:);
                Lp_move2 = likMapCell(localMap{ind(1),ind(2)}{1},W(:,mea_after_move),model);
            end
        end
    else %selected measurement from a new cell
        if isempty(mea_after_move)
            Lp_move2 = 0;
        else
            if length(mea_after_move) == 1
                Lp_move2 = Lc(mea_after_move);
            else
                in_gate = sum(gating_matrix_u(mea_after_move,:),1)>=1;
                if any(in_gate)
                    mstore = ismember(1:m,mea_after_move);
                    idx = any(all(bsxfun(@eq,reshape(mstore.',1,m,[]),Mstore),2),3);
                    if ~any(idx)
                        Lp_move2 = likPPP(W(:,mea_after_move),lik_matrix_u(mea_after_move,:),PPP,in_gate,model,filter);
                        Mstore = [Mstore;mstore];
                        Lstore = [Lstore;Lp_move2];
                    else
                        Lp_move2 = Lstore(idx);
                    end
                else
                    Lp_move2 = -inf;
                end
            end
        end
    end
    for i = 1:N
        if i==cell_idx
            Wp(i) = 0;
        else
            if i <= n_tt+n_cells
                if i<= n_tt
                    if gating_matrix_dj(mea_selected,i)
                        Lp_move1(i) = likBern(tracks(i).Bern,W(:,[meas_cell{i} mea_selected]),lik_matrix_dj([meas_cell{i} mea_selected],i),...
                            model,filter);
                    end
                else
                    if gating_matrix_m(mea_selected,i-n_tt)
                        ind = measuredCells(i-n_tt,:);
                        Lp_move1(i) = likMapCell(localMap{ind(1),ind(2)}{1},W(:,[meas_cell{i} mea_selected]),model);
                    end
                end
            else
                in_gate = sum(gating_matrix_u([meas_cell{i} mea_selected],:),1)>=1;
                if any(in_gate)
                    mstore = ismember(1:m,[meas_cell{i} mea_selected]);
                    idx = any(all(bsxfun(@eq,reshape(mstore.',1,m,[]),Mstore),2),3);
                    if ~any(idx)
                        Lp_move1(i) = likPPP(W(:,[meas_cell{i} mea_selected]),lik_matrix_u([meas_cell{i} mea_selected],:),PPP,in_gate,model,filter);
                        Mstore = [Mstore;mstore];
                        Lstore = [Lstore;Lp_move1(i)];
                    else
                        Lp_move1(i) = Lstore(idx);
                    end
                end
            end
            Wp(i) = Lp_move1(i) + Lp_move2 - Lp(cell_idx) - Lp(i);
        end
    end
    
    %selected measurement moved to a new cell
    if cell_idx <=n_tt+n_cells
        Wp(N+1) = Lp_move2 + Lc(mea_selected) - Lp(cell_idx);
    elseif cell_idx > n_tt+n_cells && ~isempty(mea_after_move)
        Wp(N+1) = Lp_move2 + Lc(mea_selected) - Lp(cell_idx);
    end
    
    if length(meas_cell{cell_idx}) > 1
        %selected cell merged with an existing cell
        Lp_merge = -inf(N,1);
        if cell_idx <=n_tt+n_cells %selected cell is a target cell
            %selected cell becomes a new cell
            in_gate = sum(gating_matrix_u(meas_cell{cell_idx},:),1)>=1;
            if any(in_gate)
                mstore = ismember(1:m,meas_cell{cell_idx});
                idx = any(all(bsxfun(@eq,reshape(mstore.',1,m,[]),Mstore),2),3);
                if ~any(idx)
                    Lp_new = likPPP(W(:,meas_cell{cell_idx}),lik_matrix_u(meas_cell{cell_idx},:),PPP,in_gate,model,filter);
                    Mstore = [Mstore;mstore];
                    Lstore = [Lstore;Lp_new];
                else
                    Lp_new = Lstore(idx);
                end
            else
                Lp_new = -inf;
            end
            Wp(N+2) = Lmiss(cell_idx) + Lp_new - Lp(cell_idx);
            
            for i = 1:n_tt %merged with a target cell
                if i~=cell_idx
                    if any(gating_matrix_dj([meas_cell{cell_idx}],i))
                        Lp_merge(i) = likBern(tracks(i).Bern,W(:,[meas_cell{cell_idx} meas_cell{i}]),lik_matrix_dj([meas_cell{cell_idx} meas_cell{i}],i),...
                            model,filter);
                    end
                    Wp(i+N+2) = Lp_merge(i) + Lmiss(cell_idx) - Lp(cell_idx) - Lp(i);
                end
            end
            for i = n_tt+1:n_cells
                if i~=cell_idx
                    if any(gating_matrix_m([meas_cell{cell_idx}],i-n_tt))
                        ind = measuredCells(i-n_tt,:);
                        Lp_merge(i) = likMapCell(localMap{ind(1),ind(2)}{1},W(:,[meas_cell{cell_idx} meas_cell{i}]),model);
                    end
                    Wp(i+N+2) = Lp_merge(i) + Lmiss(cell_idx) - Lp(cell_idx) - Lp(i);
                end
            end

            for i = n_tt+n_cells+1:N %merged with a new cell
                in_gate = sum(gating_matrix_u([meas_cell{cell_idx} meas_cell{i}],:),1)>=1;
                if any(in_gate)
                    mstore = ismember(1:m,[meas_cell{cell_idx} meas_cell{i}]);
                    idx = any(all(bsxfun(@eq,reshape(mstore.',1,m,[]),Mstore),2),3);
                    if ~any(idx)
                        Lp_merge(i) = likPPP(W(:,[meas_cell{cell_idx} meas_cell{i}]),lik_matrix_u([meas_cell{cell_idx} meas_cell{i}],:),PPP,in_gate,model,filter);
                        Mstore = [Mstore;mstore];
                        Lstore = [Lstore;Lp_merge(i)];
                    else
                        Lp_merge(i) = Lstore(idx);
                    end
                end
                Wp(i+N+2) = Lp_merge(i) + Lmiss(cell_idx) - Lp(cell_idx) - Lp(i);
            end
        else %selected cell is a new cell
            
            for i = 1:n_tt %merged with a target cell
                if any(gating_matrix_dj([meas_cell{cell_idx}],i))
                    Lp_merge(i) = likBern(tracks(i).Bern,W(:,[meas_cell{cell_idx} meas_cell{i}]),lik_matrix_dj([meas_cell{cell_idx} meas_cell{i}],i),...
                        model,filter);
                end
                Wp(i+N+2) = Lp_merge(i) - Lp(cell_idx) - Lp(i);
            end

            for i = n_tt+1:n_cells
                if any(gating_matrix_m([meas_cell{cell_idx}],i-n_tt))
                    ind = measuredCells(i-n_tt,:);
                    Lp_merge(i) = likMapCell(localMap{ind(1),ind(2)}{1},W(:,[meas_cell{cell_idx} meas_cell{i}]),model);
                end
                Wp(i+N+2) = Lp_merge(i) - Lp(cell_idx) - Lp(i);
            end
            
            for i = n_tt+n_cells+1:N %merged with a new cell
                if i~=cell_idx
                    in_gate = sum(gating_matrix_u([meas_cell{cell_idx} meas_cell{i}],:),1)>=1;
                    if any(in_gate)
                        mstore = ismember(1:m,[meas_cell{cell_idx} meas_cell{i}]);
                        idx = any(all(bsxfun(@eq,reshape(mstore.',1,m,[]),Mstore),2),3);
                        if ~any(idx)
                            Lp_merge(i) = likPPP(W(:,[meas_cell{cell_idx} meas_cell{i}]),lik_matrix_u([meas_cell{cell_idx} meas_cell{i}],:),PPP,in_gate,model,filter);
                            Mstore = [Mstore;mstore];
                            Lstore = [Lstore;Lp_merge(i)];
                        else
                            Lp_merge(i) = Lstore(idx);
                        end
                    end
                    Wp(i+N+2) = Lp_merge(i) - Lp(cell_idx) - Lp(i);
                end
            end
        end
        
        %selected cell splited into two parts
        if length(meas_cell{cell_idx}) > 2
            [IDX,~] = kmeanspp(W(:,meas_cell{cell_idx}),2);
            Lp_split1 = -inf(2,1);
            Lp_split2 = -inf(2,1);
            if cell_idx <=n_tt+n_cells %selected cell is a target cell
                if length(meas_cell{cell_idx}(IDX==2)) > 1
                    if cell_idx <= n_tt
                        Lp_split1(1) = likBern(tracks(cell_idx).Bern,W(:,meas_cell{cell_idx}(IDX==1)),lik_matrix_dj(meas_cell{cell_idx}(IDX==1),cell_idx),model,filter);
                    else
                        ind = measuredCells(cell_idx-n_tt,:);
                        Lp_split1(1) = likMapCell(localMap{ind(1),ind(2)}{1},W(:,meas_cell{cell_idx}(IDX==1)),model);
                    end
                    in_gate = sum(gating_matrix_u(meas_cell{cell_idx}(IDX==2),:),1)>=1;
                    if any(in_gate)
                        mstore = ismember(1:m,meas_cell{cell_idx}(IDX==2));
                        idx = any(all(bsxfun(@eq,reshape(mstore.',1,m,[]),Mstore),2),3);
                        if ~any(idx)
                            Lp_split2(1) = likPPP(W(:,meas_cell{cell_idx}(IDX==2)),lik_matrix_u(meas_cell{cell_idx}(IDX==2),:),PPP,in_gate,model,filter);
                            Mstore = [Mstore;mstore];
                            Lstore = [Lstore;Lp_split2(1)];
                        else
                            Lp_split2(1) = Lstore(idx);
                        end
                    end
                end
                Wp(2*N+3) = Lp_split1(1) + Lp_split2(1) - Lp(cell_idx);
                
                if length(meas_cell{cell_idx}(IDX==1)) > 1
                    if  cell_idx <= n_tt
                        Lp_split1(2) = likBern(tracks(cell_idx).Bern,W(:,meas_cell{cell_idx}(IDX==2)),lik_matrix_dj(meas_cell{cell_idx}(IDX==2),cell_idx),model,filter);
                    else
                        ind = measuredCells(cell_idx-n_tt,:);
                        Lp_split1(2) = likMapCell(localMap{ind(1),ind(2)}{1},W(:,meas_cell{cell_idx}(IDX==2)),model);
                    end
                    in_gate = sum(gating_matrix_u(meas_cell{cell_idx}(IDX==1),:),1)>=1;
                    if any(in_gate)
                        mstore = ismember(1:m,meas_cell{cell_idx}(IDX==1));
                        idx = any(all(bsxfun(@eq,reshape(mstore.',1,m,[]),Mstore),2),3);
                        if ~any(idx)
                            Lp_split2(2) = likPPP(W(:,meas_cell{cell_idx}(IDX==1)),lik_matrix_u(meas_cell{cell_idx}(IDX==1),:),PPP,in_gate,model,filter);
                            Mstore = [Mstore;mstore];
                            Lstore = [Lstore;Lp_split2(2)];
                        else
                            Lp_split2(2) = Lstore(idx);
                        end
                    end
                end
                Wp(2*N+4) = Lp_split1(2) + Lp_split2(2) - Lp(cell_idx);
                
            else %selected cell is a new cell
                if length(meas_cell{cell_idx}(IDX==1)) > 1
                    in_gate = sum(gating_matrix_u(meas_cell{cell_idx}(IDX==1),:),1)>=1;
                    if any(in_gate)
                        mstore = ismember(1:m,meas_cell{cell_idx}(IDX==1));
                        idx = any(all(bsxfun(@eq,reshape(mstore.',1,m,[]),Mstore),2),3);
                        if ~any(idx)
                            Lp_split1(1) = likPPP(W(:,meas_cell{cell_idx}(IDX==1)),lik_matrix_u(meas_cell{cell_idx}(IDX==1),:),PPP,in_gate,model,filter);
                            Mstore = [Mstore;mstore];
                            Lstore = [Lstore;Lp_split1(1)];
                        else
                            Lp_split1(1) = Lstore(idx);
                        end
                    end
                end
                if length(meas_cell{cell_idx}(IDX==2)) > 1
                    in_gate = sum(gating_matrix_u(meas_cell{cell_idx}(IDX==2),:),1)>=1;
                    if any(in_gate)
                        mstore = ismember(1:m,meas_cell{cell_idx}(IDX==2));
                        idx = any(all(bsxfun(@eq,reshape(mstore.',1,m,[]),Mstore),2),3);
                        if ~any(idx)
                            Lp_split2(1) = likPPP(W(:,meas_cell{cell_idx}(IDX==2)),lik_matrix_u(meas_cell{cell_idx}(IDX==2),:),PPP,in_gate,model,filter);
                            Mstore = [Mstore;mstore];
                            Lstore = [Lstore;Lp_split2(1)];
                        else
                            Lp_split2(1) = Lstore(idx);
                        end
                    end
                end
                Wp(2*N+3) = Lp_split1(1) + Lp_split2(1) - Lp(cell_idx);
            end
        end
    end
    
    [Wp,~] = normalizeLogWeights(Wp);
    I = find(cumsum(exp(Wp))>rand(),1); %sampling
    if I==cell_idx
        Lp = Lp_original;
        meas_cell = meas_cell_original;
        num_repetition = num_repetition + 1;
        if (num_repetition > max_repetition)
            break;
        end
    elseif I <= N %move selected measurement to an existing cell
        meas_cell{I} = [meas_cell{I} mea_selected];
        Lp(I) = Lp_move1(I);
        meas_cell{cell_idx} = mea_after_move;
        if isempty(meas_cell{cell_idx})
            if cell_idx <=n_tt+n_cells
                Lp(cell_idx) = Lmiss(cell_idx);
            else
                meas_cell(cell_idx) = [];
                Lp(cell_idx) = [];
            end
        else
            Lp(cell_idx) = Lp_move2;
        end
    elseif I==N+1 %move selected measurement to a new cell
        meas_cell{N+1,1} = mea_selected;
        Lp(N+1,1) = Lc(mea_selected);
        meas_cell{cell_idx} = mea_after_move;
        Lp(cell_idx) = Lp_move2;
    elseif I==N+2 %selected target cell becomes a new cell
        meas_cell{N+1,1} = meas_cell{cell_idx};
        Lp(N+1,1) = Lp_new;
        meas_cell{cell_idx} = [];
        Lp(cell_idx) = Lmiss(cell_idx);
    elseif I<=2*N+2 %merge selected cell with an existing cell
        meas_cell{I-N-2} = [meas_cell{cell_idx} meas_cell{I-N-2}];
        Lp(I-N-2) = Lp_merge(I-N-2);
        if cell_idx <= n_tt+n_cells
            meas_cell{cell_idx} = [];
            Lp(cell_idx) = Lmiss(cell_idx);
        else
            meas_cell(cell_idx) = [];
            Lp(cell_idx) = [];
        end
    elseif I==2*N+3 %split selected cell
        meas_cell{N+1,1} = meas_cell{cell_idx}(IDX==2);
        Lp(N+1,1) = Lp_split2(1);
        meas_cell{cell_idx} = meas_cell{cell_idx}(IDX==1);
        Lp(cell_idx) = Lp_split1(1);
    elseif I==2*N+4 %split selected cell
        meas_cell{N+1,1} = meas_cell{cell_idx}(IDX==1);
        Lp(N+1,1) = Lp_split2(2);
        meas_cell{cell_idx} = meas_cell{cell_idx}(IDX==2);
        Lp(cell_idx) = Lp_split1(2);
    end
    
    if I~=cell_idx
        num_repetition = 0;
    end
    
    partitions{t,1} = cellfun(@(x) sort(x),meas_cell,'Un',0);
    lik(t) = sum(Lp);
end

idx_empty = cellfun('isempty',partitions);
lik = lik(~idx_empty);
partitions = partitions(~idx_empty);

[lik,idx] = unique(lik);
partitions = partitions(idx); %Each partition contains a number of measurement cells

end

