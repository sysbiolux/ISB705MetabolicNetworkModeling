clear 
%% Pathway analysis
sampleModel=1; consensusModel=1;
load ('Reconstruction','ModelMatrix')

if consensusModel==1  && sampleModel==1
    name ={'consensus', 'sample'};
elseif  consensusModel==1  && sampleModel==0
    name ={'consensus'};
elseif  consensusModel==0 && sampleModel==1
    name ={'sample'};
end

if sum(ismember(name,'consensus'))
    if ~exist('ModelMatrix','var')
        sampleModelsMatrix(:,1)=ismember(ConsistentModel.rxns, model1);
        sampleModelsMatrix(:,2)=ismember(ConsistentModel.rxns, model2);% adapt here in function of the number of your models
        condition={'condition1', 'condition2'};
    end
end

ConsistentModel=ModelMatrix.consistentModel;
uniqueConditions=ModelMatrix.conditions;
[pathways, ~, ub] = unique(ConsistentModel.subSystems);
path_counts = histc(ub, 1:length(pathways));
clear ub
Pathwayskeep=struct();
Tcons = table(pathways, path_counts);
for i=1:numel(name)



    if  sum(ismember(name(i),'consensus'))
        BinaryMatrix=ModelMatrix.consensus;
        col=uniqueConditions;
    else
        BinaryMatrix=ModelMatrix.sample;
        col=ModelMatrix.colnames;
    end
    Pathways_num=zeros(size(Tcons,1),numel(col));

    %% Pathway information for the consensus models
    for counter=1:numel(col)
        [pathways, ~, ub] = unique(ConsistentModel.subSystems(find(BinaryMatrix(:,counter))));
        path_counts = histc(ub, 1:length(pathways));
        T = table(pathways, path_counts);
        [~, ia, ib] = intersect(Tcons.pathways, T.pathways);
        Pathways_num(ia, counter) = (T.path_counts(ib)) ;
    end
    Pathways=[Tcons,array2table(Pathways_num)];
    PathwayActivity = Pathways;

    Pathways.Properties.VariableNames(counter+2) = strcat(col(counter),'_',name(i));
    for counter2=3:size(PathwayActivity,2)
    PathwayActivity(:,counter2) = array2table(table2array(PathwayActivity(:,counter2))./table2array(PathwayActivity(:,2)));
    end
    Pathwayskeep(1).(name{i})=Pathways;
    Pathwayskeep(2).(name{i})=PathwayActivity;
    clear T pathways
end
clear i col Pathways Pathways_num ia ib counter counter2 BinaryMatrix
