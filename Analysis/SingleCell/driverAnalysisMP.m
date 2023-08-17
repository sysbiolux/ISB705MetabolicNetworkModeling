% driver analysis, TS 07/23
%% scFASTCORMICS HPC results Maria 10.07.23
clearvars -except solverOK
clc, close all
close force all

%% Load model/define the clusters. conditions and the medium

% load('scforClaudia_no_medium.mat')
load('Claudia_medium_cons_250723.mat')
% load('Claudia_model_wMediumCons_270723.mat')
clearvars -except solverOK Results_keep
load('simpleRev_recon3.mat')
Recon=model;
load('Table_medium.mat')
T.medium_rxns_keep(cellfun('isempty',T.medium_rxns_keep))=cellstr('nd');
T=[T; {'Ex_dttp[e]','dttp[e]',1}];
ClusterNames={'c1';'c2'; 'c3'; 'c4';'c5'};
names_col={'Ctrl30';'Ctrl60'; 'PD30'; 'PD60';'GC30'; 'GC60'};
changeCobraSolver('ibm_cplex')
clear ans

%% model sizes
Modelsizes=zeros(size(Results_keep,1),1);
for counter=1:numel(Results_keep)
    Modelsizes(counter)=numel(Results_keep(counter).multi_cell_population.rxns);
end
disp('Model sizes:')
disp(names_col')
disp(Modelsizes)
clear counter

%% model sizes per cluster
clusterSize=zeros(size(Results_keep,2),numel(ClusterNames));
for counter=1:size(clusterSize,1)
    %     model=Results_keep(counter).multi_cell_population_model;
    model=Results_keep(counter).multi_cell_population;
    for counter2=1:size(clusterSize,2)
        idx = strfind(model.rxns, ['_' num2str(counter2)]);
        idx = find(not(cellfun('isempty', idx)));
        clusterSize(counter,counter2)=numel(idx);
    end
end
disp('Model sizes per clusters (cols):')
clusterSize=array2table(clusterSize);
clusterSize.Properties.VariableNames={'c1','c2','c3','c45','c6'};
clear counter counter2 idx ClusterNames
clusterSize.Row=names_col;
disp(clusterSize)


%% Jaccard plot Whole Models
J=nan(numel(Results_keep));
for counter=1:numel(Results_keep)
    for counter2=1:numel(Results_keep)
        A1=find(Results_keep(counter).A);
        A2=find(Results_keep(counter2).A);
        J(counter,counter2)=numel(intersect(A1,A2))/numel(union(A1,A2));
    end
end
disp('Jaccard similarity Whole Models:')

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
    'Colormap', altcolor);
addTitle(cgo_J,{'Jaccard similarity Whole Models'});
plot(cgo_J);
clear counter counter2 J A1 A2

%% Jaccard plot per Cluster over all models

%%% Why are reactions in each cluster in generic model slightly different?

for counter=1:numel(clusterSize.Properties.VariableNames) %per cluster
    J=nan(numel(clusterSize.Properties.VariableNames));
    % identify from generic model rxns in specific cluster
    model=Results_keep(1).ExpandedInputModel;
    idx = strfind(model.rxns, ['_' num2str(counter)]);
    idx = find(not(cellfun('isempty', idx)));

    for counter2=1:numel(Results_keep) %per model
        for counter3=1:numel(Results_keep) %per model
            A1=find(Results_keep(counter2).A);
            A1=intersect(A1,idx);
            A2=find(Results_keep(counter3).A);
            A2=intersect(A2,idx);
            J(counter2,counter3)=numel(intersect(A1,A2))/numel(union(A1,A2));
        end
    end
    disp(['Jaccard similarity Whole Cluster:' num2str(counter)])

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
        'Colormap', altcolor);
    addTitle(cgo_J,{['Jaccard similarity Cluster: ' num2str(counter)]});
    plot(cgo_J);
end
clear A1 A2 idx counter counter2 counter3 J

%% Jaccard plot per all Clusters within one model
for counter=1:numel(Results_keep) %per model
    J=nan(numel(clusterSize.Properties.VariableNames));
    %     model=Results_keep(counter).multi_cell_population_model;
    model=Results_keep(counter).multi_cell_population;

    for counter2=1:numel(clusterSize.Properties.VariableNames) %per cluster
        idx1 = strfind(model.rxns, ['_' num2str(counter2)]);
        idx1 = find(not(cellfun('isempty', idx1)));
        A1 = cellfun(@(S) S(1:end-2), model.rxns(idx1), 'Uniform', 0);
        for counter3=1:numel(clusterSize.Properties.VariableNames) %per cluster
            idx2 = strfind(model.rxns, ['_' num2str(counter3)]);
            idx2 = find(not(cellfun('isempty', idx2)));
            A2 = cellfun(@(S) S(1:end-2), model.rxns(idx2), 'Uniform', 0);
            J(counter2,counter3)=numel(intersect(A1,A2))/numel(union(A1,A2));
        end
    end
    disp(['Jaccard similarity Model:' num2str(counter)])

    altcolor = [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;...
        255 0 0; 204 0 0; 152 0 0; 102 0 0; 51 0 0]/255;
    cgo_J = clustergram(J,...
        'RowLabels', clusterSize.Properties.VariableNames,...
        'ColumnLabels', clusterSize.Properties.VariableNames,...
        'ColumnLabelsRotate',340, ...
        'Cluster', 'all', ...
        'Annotate', 'true',...
        'symmetric','False',...
        'AnnotColor','k',...
        'Colormap', altcolor);
    addTitle(cgo_J,{['Jaccard similarity Model: ' num2str(counter)]});
    plot(cgo_J);
end

%% Pathway anaysis model centred
multiCellModel=Results_keep(1).multi_cell_population;
ExpandedInputModel=Results_keep(1).ExpandedInputModel;
clear counter counter2 counter3 idx1 idx2 J A1 A2
close all

Pathways_keep=struct();
for i=1:numel(names_col)
    multiCellModel=Results_keep(i).multi_cell_population;
    [Pathways]=PathwayAnalysisMultiCellModel(multiCellModel,ExpandedInputModel,  clusterSize.Properties.VariableNames);

    Pathways_keep(i).Analysis=Pathways;
    Pathways_keep(i).multiCellModel=multiCellModel;
    %Making a plot that show interesting pathways differences
    %between the different clusters (Tumor)

    %Finding the size of the top 30 most differentially present pathways
    Diff=zeros(size(Pathways,1),1);
    %Finding top thirty pathways differential
    for ii = 1:height(Pathways)
        maxValue = max(Pathways{ii,:});
        minValue = min(Pathways{ii,:});
        Diff(ii,1) = (maxValue - minValue);
    end
    Pathways2=[Pathways, table(Diff)];
    index=find(~isnan(Pathways2.Diff));

    ClusterRatioShort=Pathways2(index(1:30),:);

    cgo_Path = clustergram(table2array(ClusterRatioShort(:,1:end-1)),...
        'RowLabels',ClusterRatioShort.Row ,...
        'ColumnLabels',ClusterRatioShort.Properties.VariableNames(1:end-1),...
        'ColumnLabelsRotate',340, ...
        'Cluster', 'all',...
        'LabelsWithMarkers',true,...
        'Colormap',redbluecmap,...
        'DisplayRange',1000,...
        'DisplayRatio', [0.2 0.2]);

    addTitle(cgo_Path, cellstr(strcat('Pathway Analysis',' ', names_col{i})));

    graph = plot(cgo_Path);
    set(graph, 'Clim', [0,1])
    colormap(graph,"copper")

    figureHandle = gcf;
    %# make all text in the figure to wanted size and bold
    set(findall(figureHandle,'type','text'),'fontSize',16,'fontWeight','bold')
    set(findall(figureHandle, 'type','axes','tag','HeatMapAxes'),'fontsize',16)

end


%% Pathway anaysis cluster centred
[ClusterRatioKeep]=PathwayAnalysisMultiCellModel4Cluster(Results_keep,ExpandedInputModel, clusterSize.Properties.VariableNames, names_col);
for ii=1:numel(ClusterRatioKeep)
    % Finding the size of the top 30 most differentially present pathways
    ClusterRatio=ClusterRatioKeep(ii).cluster;
    Diff=zeros(size(ClusterRatio,1),1);

    for i=1:size(ClusterRatio,1)
        %Finding top thirty pathways differential
        maxValue = max(ClusterRatio{i,:});
        minValue = min(ClusterRatio{i,:});
        Diff(i,1) = (maxValue - minValue);
    end
    ClusterRatio2=[ClusterRatio, table(Diff)];
    index=find(~isnan(ClusterRatio2.Diff));

    ClusterRatioShort=ClusterRatio2(index(1:30),:);


    cgo_Path = clustergram(table2array(ClusterRatioShort(:,1:end-1)),...
        'RowLabels',ClusterRatioShort.Row ,...
        'ColumnLabels',ClusterRatioShort.Properties.VariableNames(1:end-1),...
        'ColumnLabelsRotate',340, ...
        'Cluster', 'all',...
        'LabelsWithMarkers',true,...
        'Colormap',redbluecmap,...
        'DisplayRange',1000,...
        'DisplayRatio', [0.2 0.2]);

    addTitle(cgo_Path, cellstr(strcat('Pathway Analysis',' ', clusterSize.Properties.VariableNames{ii})));

    graph = plot(cgo_Path);
    set(graph, 'Clim', [0,1])
    colormap(graph,"copper")

    figureHandle = gcf;
    %# make all text in the figure to wanted size and bold
    set(findall(figureHandle,'type','text'),'fontSize',16,'fontWeight','bold')
    set(findall(figureHandle, 'type','axes','tag','HeatMapAxes'),'fontsize',16)
end
close all


%% run this for all 6 models ... (change index in script!)
% FBA
sortFlag=1; %1:max-min; 2:max; 3:min

[AnalysisKeep,ExchangesKeep]=ExchangesPrediction(Results_keep, Recon, T, clusterSize);

for i=1:numel(ExchangesKeep)
    ClusterExchanges=ExchangesKeep(i).Exchanges;
    sum_of_matrix_rows = sum(table2array(ClusterExchanges),2)~=0;
    ClusterExchangesShort=ClusterExchanges(sum_of_matrix_rows,:);



    cgo_Path = clustergram(table2array(ClusterExchangesShort),...
        'RowLabels',ClusterExchangesShort.Row ,...
        'ColumnLabels', ClusterExchangesShort.Properties.VariableNames,...
        'ColumnLabelsRotate',340, ...
        'Cluster', 'all',...
        'LabelsWithMarkers',true,...
        'Colormap',redbluecmap,...
        'DisplayRange',10,...
        'DisplayRatio', [0.2 0.2]);

    addTitle(cgo_Path, cellstr(strcat('FBA of exchanged metabolites','Model:', names_col(i), ', sortflag: ',num2str(sortFlag))))

    plot(cgo_Path);

    figureHandle = gcf;
    %# make all text in the figure to size 14 and bold
    set(findall(figureHandle,'type','text'),'fontSize',16,'fontWeight','bold')
    set(findall(figureHandle, 'type','axes','tag','HeatMapAxes'),'fontsize',16)

    % red (positive values) indicates the metabolite leaves the corresponding cluster
    % blue (negative values) indicates the metabolite enters the corresponding cluster
    fileName=['data' num2str(i) '.xlsx'];
    if exist(fileName, 'file') == 2

        delete(fileName)
    end
    writetable(ExchangesKeep(i).Exchanges,fileName,'WriteRowNames',1)
end

%save(['FBAmodel' num2str(modelNr)],'FBAweighted','multiCellModel')
[ComparisonFlux2Models,res4t, FluxPerCluster,FluxPerPathway, FluxPerPathwayPerCluster]=compareFBA(AnalysisKeep, Results_keep);
k=0;
for i=1:numel(FluxPerCluster)
    disp('reactions carrying a flux per cluster')
    disp(FluxPerCluster(i).FluxPerCluster)
end
for i=1:numel(FluxPerPathway)

    disp('reactions carrying a flux per pathways')
    disp(FluxPerPathway)
end
for counterA=1:numel(Results_keep)
    for counterB=1:numel(Results_keep)
         if counterA<counterB 

    k=k+1;

    disp('relative flux difference per subSystem and cluster')
    res4= res4t(k).rest;
    res4m=table2array(res4);
    temp=sort(max(res4m,[],2),'descend');
    show=find(max(res4m,[],2)>=temp(10)); %show top10
    resTop10=res4(show,:);


    fileName=['compareFBA_' num2str(counterA) '_vs' num2str(counterB) '.xlsx'];
    if exist(fileName, 'file') == 2

        delete(fileName)
    end
    writetable(res4,fileName,'WriteRowNames',1)
    end
    end
end
