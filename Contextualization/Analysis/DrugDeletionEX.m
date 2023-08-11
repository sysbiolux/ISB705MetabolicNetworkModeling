clear 
load GeneDrugRelations.mat
DrugList = unique(GeneDrugRelations.DrugName)
sampleModel=1; consensusModel=1;
load ('Reconstruction','ModelMatrix')
changeCobraSolver('ibm_cplex')
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

objectives={'DM_atp_c_', 'biomass_reaction'};
DrugResults=struct();
for counter2=1:numel(name)

    if  sum(ismember(name(counter2),'consensus'))
        models_keep=ModelMatrix.consensus;
    else
        models_keep=ModelMatrix.sample;
    end
    for i=1:size(models_keep,2) %for each model

        for counter=1:numel(objectives)

            ind = find(~cellfun(@isempty, regexp(ConsistentModel.rxns, objectives (counter) )));
            model_out = removeRxns(ConsistentModel,ConsistentModel.rxns(setdiff(1:numel(ConsistentModel.rxns),find(models_keep(:,i))))); % create model based on active reactions
            if sum(ismember(model_out.rxns,objectives(counter)))
                model_out = changeObjective(model_out,ConsistentModel.rxns(ind));


               [grRatio, grRateKO, grRateWT] = DrugDeletion(model_out,'FBA',DrugList);

            DrugResults(i).(name{counter2}). (objectives{counter}).gratio=grRatio;
            DrugResults(i). (name{counter2}). (objectives{counter}).grRateKO=grRateKO;
            DrugResults(i).(name{counter2}). (objectives{counter}).grRateWT=grRateWT;
            else
                disp( [objectives(counter),'is not present in model', num2str(i)] )
            end
        end
    end
end

save DrugDeletion