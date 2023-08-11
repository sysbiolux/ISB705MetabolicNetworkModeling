% sampleModelsMatrix= binary matrix with n rows corresponding to the number of reactions and m columns corresponding
%to the number of models

altcolor= [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;...
    255 0 0; 204 0 0; 152 0 0; 102 0 0;  51 0 0]/255; %shorter 10% = 1 bar
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
    sampleModelsMatrix(:,1)=ismember(ConsistentModel.rxns, model1)
    sampleModelsMatrix(:,2)=ismember(ConsistentModel.rxns, model2)% adapt here in function of the number of your models
    condition={'condition1', 'condition2'}
end
end

for i=1:numel(name)
    if  sum(ismember(name(i),'consensus'))
        col=ModelMatrix.conditions;
        modelkeep=ModelMatrix.consensus;
    else
        col=ModelMatrix.colnames
        modelkeep=ModelMatrix.sample;
    end
% similarity between two models can be assessed via the Jaccard similarity index
J = squareform(pdist(modelkeep','jaccard'));
cgo_J = clustergram(1-J,...
    'RowLabels', col,...
    'ColumnLabels', [],...
    'ColumnLabelsRotate',270, ...
    'Cluster', 'all', ...
    'symmetric','False',...
    'Colormap', altcolor);
addTitle(cgo_J,{'Model similarity based on Jaccard distance','sampleModelsMatrix'})
end