clear
load GeneDrugRelations.mat
DrugList = unique(GeneDrugRelations.DrugName);
sampleModel=1; consensusModel=1;
changeCobraSolver('ibm_cplex')

if consensusModel==1  && sampleModel==1
    annotation ={'consensus', 'sample'};
elseif  consensusModel==1  && sampleModel==0
    annotation ={'consensus'};
elseif  consensusModel==0 && sampleModel==1
    annotation ={'sample'};
end


objectives={'DM_atp_c_', 'biomass_reaction'};
DrugResults=struct();

for counter2=1:numel(annotation)
    
    if  sum(ismember(annotation(counter2),'consensus'))
        load Reconstruction uniqueconditions
        col=uniqueconditions;
        name=strcat('ConsensusModel_',uniqueconditions{counter2});
        load (name, 'ContextModel');
        if ischar(name)
            name=cellstr(name);
        end
        
    else
        load Reconstruction colnames
        
        name=strcat('SampleModel_',colnames{counter2});
        
        load (name,'ContextModel');
        col=colnames;
        if ischar(name)
            name=cellstr(name);
        end
        
    end
    
    for i=1:2%numel(col) %for each model
        
        for counter=1:numel(objectives)
            
            if sum(ismember(ContextModel.rxns,objectives(counter)))
                ContextModel.c=zeros(numel(ContextModel.rxns),1);
                ContextModel.c(ismember(ContextModel.rxns,objectives(counter)))=1;
                
                [grRatio, grRateKO, grRateWT] = DrugDeletion( ContextModel,'FBA',DrugList);
                
                DrugResults(i).(annotation{counter}). (objectives{counter}).gratio=grRatio;
                DrugResults(i). (annotation{counter}). (objectives{counter}).grRateKO=grRateKO;
                DrugResults(i).(annotation{counter}). (objectives{counter}).grRateWT=grRateWT;
            else
                disp( [objectives(counter),'is not present in model', num2str(i)] )
            end
        end
    end
end

save DrugDeletion