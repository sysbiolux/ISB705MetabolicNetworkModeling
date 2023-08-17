function[ClusterRatioKeep]=PathwayAnalysisMultiCellModel4Cluster(Results_keep,ExpandedInputModel, cluster, ModelNames)
%Pathways Analysis Code for scFASTCORMICS output (models with
%strictly less than 100 clusters)

%(c) Leo Eusebi, Luc Wagener 2023 et al -University of Luxembourg
%%

%% All (specific) models for a specific cluster
%creating a string of the existing unique pathways
Pathways = unique(string(ExpandedInputModel.subSystems));
Pathways(ismember(Pathways,"")) = [];
%creating a matrix with pathways as rows and clusters as columns, the
%values shown in that matrix will be the quantity of each pathways
%depending on each clusters
ClusterRatioKeep=struct();
Cluster = cell(numel(ExpandedInputModel.rxns),1);
LastTwo = cellfun(@(S) S(end-1), ExpandedInputModel.rxns, 'Uniform', 0);
match=contains(LastTwo,'_');
Cluster(match) = cellfun(@(S) S(end), ExpandedInputModel.rxns(match), 'Uniform', 0);
match=~contains(LastTwo,'_') & contains(LastTwo,'e');
Cluster(match)=cellstr('e');
if cellfun('isempty',Cluster)
    save
sss
end

for counter=1:numel(cluster)
    ClusterMat=zeros(numel(Pathways),numel(ModelNames));
ClusterMatExpanded=zeros(numel(Pathways),numel(ModelNames));

for i = 1:size(ModelNames,1)
          ind=ismember(ExpandedInputModel.rxns,Results_keep(i).multi_cell_population.rxns);
%         IdxCluster = find(contains(LastTwo,strcat('_',num2str(i))));
        IdxCluster = find(contains(Cluster(ind),num2str(counter)));   
   
    for u = 1:length(Pathways)
        
        PathwaysInModels = ismember(string(Results_keep(i).multi_cell_population.subSystems(IdxCluster)),Pathways(u)) ;
        ClusterMat(u,i) = sum(PathwaysInModels);
        PathwaysInModelsEx = ismember(string(ExpandedInputModel.subSystems(IdxCluster)),Pathways(u)) ;
        ClusterMatExpanded(u,i) = sum(PathwaysInModelsEx);        
    end
end


ClusterMat = array2table(ClusterMat);
ClusterMat.Properties.RowNames = cellstr(Pathways) ;
ClusterMat.Properties.VariableNames = ModelNames ;

ClusterMatExpanded = array2table(ClusterMatExpanded);
ClusterMatExpanded.Properties.RowNames = cellstr(Pathways) ;
ClusterMatExpanded.Properties.VariableNames = ModelNames;

%%


%Standarization to help in the visualisation process on the plot
ClusterRatio = table2array(ClusterMat)./table2array(ClusterMatExpanded(:,1));
ClusterRatio = array2table(ClusterRatio);
ClusterRatio.Properties.RowNames = Pathways(:) ;

ClusterRatio.Properties.VariableNames = ModelNames; %PDNames ;
%%

%Making a plot that show interesting pathways differences
%between the different clusters (Tumor)

ClusterRatioKeep(counter).cluster=ClusterRatio;



end
end