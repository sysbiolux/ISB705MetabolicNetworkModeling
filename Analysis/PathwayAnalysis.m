%% Pathway analysis
clear
load ('Reconstruction','SampleMatrix','colnames','ConsistentModel','ConsensusMatrix')
altcolor= [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;...
    255 0 0; 204 0 0; 152 0 0; 102 0 0;  51 0 0]/255; %shorter 10% = 1 bar

if iscell(ConsistentModel.subSystems{1})
    ConsistentModel.subSystems=vertcat(ConsistentModel.subSystems{:})
end
Pathways = table(unique(ConsistentModel.subSystems));
[pathways, ~, ub] = unique(ConsistentModel.subSystems);
path_counts = histc(ub, 1:length(pathways));
T = table(pathways, path_counts);
[~, ia, ib] = intersect(Pathways.Var1, T.pathways);
Pathways.consistent(ia) = T.path_counts(ib);
%% Pathway information for the consensus models

[pathways, ~, ub] = unique(ConsistentModel.subSystems(find(ConsensusMatrix(:,1))));
path_counts = histc(ub, 1:length(pathways));
T = table(pathways, path_counts);
[~, ia, ib] = intersect(Pathways.Var1, T.pathways);
Pathways.Var2(ia) = T.path_counts(ib) ;
Pathways.Properties.VariableNames{3} = 'cancer_consensus';
[pathways, ~, ub] = unique(ConsistentModel.subSystems(find(ConsensusMatrix(:,2))));
path_counts = histc(ub, 1:length(pathways));
T = table(pathways, path_counts);
[~, ia, ib] = intersect(Pathways.Var1, T.pathways);
Pathways.Var2(ia) = T.path_counts(ib) ;
Pathways.Properties.VariableNames{4} = 'control_consensus';
%% pathway information for the sample-specific models

for i=1:numel(colnames)
    [pathways, ~, ub] = unique(ConsistentModel.subSystems(find(SampleMatrix(:,i))));
    path_counts = histc(ub, 1:length(pathways));
    T = table(pathways, path_counts);
    [~, ia, ib] = intersect(Pathways.Var1, T.pathways);
    Pathways.Var2(ia) = T.path_counts(ib) ;
    Pathways.Properties.VariableNames{4+i} = colnames{i};
end
%% pathway activity rates

PathwayActivity = Pathways;
for i=3:size(PathwayActivity,2)
    PathwayActivity(:,i) = array2table(table2array(PathwayActivity(:,i))./table2array(PathwayActivity(:,2)));
end
% comparison of 2 conditions
% pathways with a difference higher than 20% 
diff_idx = find(abs(table2array(PathwayActivity(:,3))- table2array(PathwayActivity(:,4))) > 0.2);
%% plotting

figure
hold on
scatter(table2array(PathwayActivity(:,3)),table2array(PathwayActivity(:,4)),'filled',...
    'MarkerFaceColor',[0.9 0.9 0.9])
scatter(table2array(PathwayActivity(diff_idx,3)),table2array(PathwayActivity(diff_idx,4)),...
    'black')
ylabel('cancer consensus model')
xlabel('control consensus model')
title('Pathway presence rate in the consensus models')
line([0 1], [0,1],'Color','k')
line([0 0.8], [0.2,1],'Color','k','LineStyle','--')
line([0.2 1], [0,0.8],'Color','k','LineStyle','--')
legend({'All pathways','>20%'},"Location","best")
text(table2array(PathwayActivity(diff_idx,3)),table2array(PathwayActivity(diff_idx,4)), PathwayActivity.Var1(diff_idx))
%% comparison of multiple samples

[~,I] = sort(sum(abs(table2array(PathwayActivity(:,3:end))-mean(table2array(PathwayActivity(:,3:end)),2)),2),'descend');
cgo = clustergram(table2array(PathwayActivity(I(1:20),3:end)),...
    'RowLabels', regexprep(PathwayActivity.Var1(I(1:20)),'metabolism',''),...
    'ColumnLabels', regexprep(PathwayActivity.Properties.VariableNames(3:end),'_TCGA.*',''),...
    'ColumnLabelsRotate',270, ...
    'Cluster', 'all', ...
    'symmetric','False',...
    'Colormap', altcolor);
addTitle(cgo,'Pathway activity for all models');
h = plot(cgo); set(h,'TickLabelInterpreter','none');
colorbar(h)
save Pathway