clearvars -except solverOK
clc, close all
% load('scforClaudia_no_medium.mat')
load('\\atlas.uni.lux\Fstc_sysbio\6-EXCHANGE\Claudia_2023\scFASTCORMICS\Claudia_medium_cons_250723.mat', 'Results_keep')


%% model sizes
res =zeros(numel(Results_keep),2;
for counter=1:numel(Results_keep)
%     model=Results_keep(counter).multi_cell_population_model;
    model=Results_keep(counter).multi_cell_population;
    res=[res; counter numel(model.rxns)];
end
disp('Model sizes:')
disp(res)

%% model sizes per cluster
res=[];
for counter=1:numel(Results_keep)
%     model=Results_keep(counter).multi_cell_population_model;
    model=Results_keep(counter).multi_cell_population;
    temp=[];
    for counter2=1:5
        idx = strfind(model.rxns, ['_' num2str(counter2)]);
        idx = find(not(cellfun('isempty', idx)));
        temp=[temp, numel(idx)];
    end
    res=[res; temp];
end
disp('Model sizes per clusters (cols):')
disp(res)
namesCluster={'c1','c2','c3','c45','c6'}
 
%% Jaccard plot Whole Models
res=nan(numel(Results_keep));
for counter=1:numel(Results_keep)
    for counter2=1:numel(Results_keep)
        A1=find(Results_keep(counter).A);
        A2=find(Results_keep(counter2).A);
        res(counter,counter2)=numel(intersect(A1,A2))/numel(union(A1,A2));
    end
end
disp('Jaccard similarity Whole Models:')
disp(res)

J=res;
names_col={'Ctrl30';'Ctrl60'; 'PD30'; 'PD60';'GC30'; 'GC60'};
altcolor = [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;...
    255 0 0; 204 0 0; 152 0 0; 102 0 0; 51 0 0]/255;
cgo_J = clustergram(J,...
    'RowLabels', names_col,...
    'ColumnLabels', names_col,...
    'ColumnLabelsRotate',340, ...
    'Cluster', 'all', ...
    'Annotate', 'true',...
    'symmetric','False',...
    'AnnotColor','k',...
    'Colormap', altcolor)
addTitle(cgo_J,{'Jaccard similarity Whole Models'});
plot(cgo_J);

%% Jaccard plot per Cluster over all models

%%% Why are reactions in each cluster in generic model slightly different?

nrClusters=5;
for counter=1:nrClusters %per cluster
    counter
    res=nan(nrClusters);
    % identify from generic model rxns in specific cluster
    model=Results_keep(1).ExpandedInputModel;
    idx = strfind(model.rxns, ['_' num2str(counter)]);
    idx = find(not(cellfun('isempty', idx)));
    numel(idx)
    
    for counter2=1:numel(Results_keep) %per model
        for counter3=1:numel(Results_keep) %per model
            %             counter3
            A1=find(Results_keep(counter2).A);
            %             numel(A1)
            A1=intersect(A1,idx);
            %             numel(A1)
            A2=find(Results_keep(counter3).A);
            A2=intersect(A2,idx);
            res(counter2,counter3)=numel(intersect(A1,A2))/numel(union(A1,A2));
        end
    end
    disp(['Jaccard similarity Whole Cluster:' num2str(counter)])
    disp(res)
    
    J=res;
    names_col={'Ctrl30';'Ctrl60'; 'PD30'; 'PD60';'GC30'; 'GC60'};
    altcolor = [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;...
        255 0 0; 204 0 0; 152 0 0; 102 0 0; 51 0 0]/255;
    cgo_J = clustergram(J,...
        'RowLabels', names_col,...
        'ColumnLabels', names_col,...
        'ColumnLabelsRotate',340, ...
        'Cluster', 'all', ...
        'Annotate', 'true',...
        'symmetric','False',...
        'AnnotColor','k',...
        'Colormap', altcolor)
    addTitle(cgo_J,{['Jaccard similarity Cluster: ' num2str(counter)]});
    plot(cgo_J);
end

%% Jaccard plot per all Clusters within one model
nrClusters=5;
for counter=1:numel(Results_keep) %per model
    res=nan(nrClusters);
%     model=Results_keep(counter).multi_cell_population_model;
    model=Results_keep(counter).multi_cell_population;
    
    for counter2=1:nrClusters %per cluster
        idx1 = strfind(model.rxns, ['_' num2str(counter2)]);
        idx1 = find(not(cellfun('isempty', idx1)));
        A1 = cellfun(@(S) S(1:end-2), model.rxns(idx1), 'Uniform', 0);
        for counter3=1:nrClusters %per cluster
            idx2 = strfind(model.rxns, ['_' num2str(counter3)]);
            idx2 = find(not(cellfun('isempty', idx2)));
            A2 = cellfun(@(S) S(1:end-2), model.rxns(idx2), 'Uniform', 0);
            res(counter2,counter3)=numel(intersect(A1,A2))/numel(union(A1,A2));
        end
    end
    disp(['Jaccard similarity Model:' num2str(counter)])
    disp(res)
    
    J=res;
    names_col={'c1';'c2'; 'c3'; 'c45';'c6'};
    altcolor = [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;...
        255 0 0; 204 0 0; 152 0 0; 102 0 0; 51 0 0]/255;
    cgo_J = clustergram(J,...
        'RowLabels', names_col,...
        'ColumnLabels', names_col,...
        'ColumnLabelsRotate',340, ...
        'Cluster', 'all', ...
        'Annotate', 'true',...
        'symmetric','False',...
        'AnnotColor','k',...
        'Colormap', altcolor)
    addTitle(cgo_J,{['Jaccard similarity Model: ' num2str(counter)]});
    plot(cgo_J);
end

%% run this for all 6 models ... (change index in script!)
PathwayAnalysis

%% run this for all 5 clusters ... (change index in script!)
PathwayAnalysis2

%% run this for all 6 models ... (change index in script!)
% FBA
FBA_mediumConc
compareFBA
