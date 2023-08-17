function[ClusterRatio]=PathwayAnalysisMultiCellModel(multiCellModel,ExpandedInputModel,  ClusterNames)
%Pathways Analysis Code for scFASTCORMICS output (models with
%strictly less than 100 clusters)

%(c) Leo Eusebi, Luc Wagener 2023 et al -University of Luxembourg
%% 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Specific model
%creating a string of the existing unique pathways
Pathways = unique(string(ExpandedInputModel.subSystems));
Pathways(ismember(Pathways,"")) = [];
%creating a matrix with pathways as rows and clusters as columns, the
%values shown in that matrix will be the quantity of each pathways
%depending on each clusters
Cluster = cell(numel(multiCellModel.rxns),1);
LastTwo = cellfun(@(S) S(end-1), multiCellModel.rxns, 'Uniform', 0);
match=contains(LastTwo,'_');
Cluster(match) = cellfun(@(S) S(end), multiCellModel.rxns(match), 'Uniform', 0);
match=~contains(LastTwo,'_') & contains(LastTwo,'e');
Cluster(match)=cellstr('e');
if cellfun('isempty',Cluster)
    save
sss
end

%%
ClusterMat=zeros(numel(Pathways),numel(unique(Cluster))-1);
ClusterMatExpanded=zeros(numel(Pathways),numel(unique(Cluster))-1);

for i = 1:numel(ClusterNames) 
    IdxCluster = ismember(Cluster,num2str(i));
   

    for u = 1:length(Pathways)
        
        PathwaysInModels = ismember(string(multiCellModel.subSystems(IdxCluster)),Pathways(u)) ;
        ClusterMat(u,i) = sum(PathwaysInModels);
        PathwaysInModelsEx = ismember(string(ExpandedInputModel.subSystems(IdxCluster)),Pathways(u)) ;
        ClusterMatExpanded(u,i) = sum(PathwaysInModelsEx);

        
    end
    
end


ClusterMat = array2table(ClusterMat);
ClusterMat.Properties.RowNames = cellstr(Pathways) ;
ClusterMat.Properties.VariableNames = ClusterNames ;

ClusterMatExpanded = array2table(ClusterMatExpanded);
ClusterMatExpanded.Properties.RowNames = cellstr(Pathways) ;
ClusterMatExpanded.Properties.VariableNames = ClusterNames;

%%


%Standarization to help in the visualisation process on the plot
ClusterRatio = table2array(ClusterMat)./table2array(ClusterMatExpanded(:,1));
ClusterRatio = array2table(ClusterRatio);
ClusterRatio.Properties.RowNames = Pathways(:) ;

ClusterRatio.Properties.VariableNames = ClusterNames; %PDNames ;
end
%%


