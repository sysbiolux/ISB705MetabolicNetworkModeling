% sampleModelsMatrix= binary matrix with n rows corresponding to the number of reactions and m columns corresponding
%to the number of models
clear
altcolor= [255 255 255;255 204 204; 255 153 153; 255 102 102; 255 51 51;...
    255 0 0; 204 0 0; 152 0 0; 102 0 0;  51 0 0]/255; %shorter 10% = 1 bar
sampleModel=1; consensusModel=1;
if consensusModel==1  && sampleModel==1
    annotation ={'consensus', 'sample'};
    load Reconstruction ConsensusMatrix SampleMatrix colnames uniqueconditions
elseif  consensusModel==1  && sampleModel==0
    annotation ={'consensus'};
    load Reconstruction ConsensusMatrix uniqueconditions
elseif  consensusModel==0 && sampleModel==1
    annotation ={'sample'};
    load Reconstruction SampleMatrix colnames
end

if sum(ismember(annotation,'consensus'))
    
    if ~exist('ConsensusMatrix','var')
        load Reconstruction uniqueconditions ConsistentModel ConsensusModel
        
        ConsensusMatrix=zeros(numel(ConsistentModel.rxns),numel(uniqueconditions));
        
        for i=1:numel(ConsensusModel)
            model=ConsensusModel(i).model;
            ConsensusMatrix(:,i)=ismember(ConsistentModel.rxns, model.rxns);
        end
    end
elseif ~exist('SampleMatrix','var')
    load Reconstruction colnames SampleModel
    
    SampleMatrix=zeros(numel(ConsistentModel.rxns),numel(colnames));
    
    for i=1:numel(SampleModel(i).model)
        model=SampleModel(i).model;
        SampleMatrix(:,i)=ismember(ConsistentModel.rxns, model);
    end
end

for i=1:numel(annotation)
    if  sum(ismember(annotation(i),'consensus'))
        col=uniqueconditions;
        modelkeep=ConsensusMatrix;
    else
        col=colnames;
        modelkeep=SampleMatrix;
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