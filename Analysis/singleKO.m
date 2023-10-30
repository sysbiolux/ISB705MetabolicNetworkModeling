clear
sampleModel=1; consensusModel=1;
load ('Reconstruction','ModelMatrix')
changeCobraSolver('ibm_cplex')

if consensusModel==1  && sampleModel==1
    annotation ={'consensus', 'sample'};
elseif  consensusModel==1  && sampleModel==0
    annotation ={'consensus'};
elseif  consensusModel==0 && sampleModel==1
    annotation ={'sample'};
end

if sum(ismember(annotation,'consensus'))
    if ~exist('ModelMatrix','var')
        sampleModelsMatrix(:,1)=ismember(ConsistentModel.rxns, model1);
        sampleModelsMatrix(:,2)=ismember(ConsistentModel.rxns, model2);% adapt here in function of the number of your models
        condition={'condition1', 'condition2'};
    end
end

objectives={'DM_atp_c_', 'biomass_reaction'};
SingleKOResults=struct();
for counter2=1:numel(annotation)

    if  sum(ismember(annotation(counter2),'consensus'))
        load Reconstruction uniqueconditions
        col=uniqueconditions;
        name=strcat('ConsensusModel_',uniqueconditions{counter2});
        load (name, 'ContextModel');
        if ischar(name)
            name=cellstr(name)
        end
        
    else
         load Reconstruction colnames

        name=strcat('SampleModel_',colnames{counter2});
        load (name, 'ContextModel');
        col=colnames;
        if ischar(name)
            name=cellstr(annotation)
        end
        
    end
    for i=1:numel(col) %for each model

        for counter=1:numel(objectives)


            if sum(ismember(ContextModel.rxns,objectives(counter)))
                 ContextModel.c=zeros(numel(ContextModel.rxns),1);
                ContextModel.c(ismember(ContextModel.rxns,objectives(counter)))=1;
                [grRatio, grRateKO, grRateWT, hasEffect, delRxns, fluxSolution, geneList] = singleGeneDeletion_rFASTCORMICS(ContextModel,'FBA',[],0,1);

                SingleKOResults(i).(name{counter2}). (objectives{counter}).gratio=grRatio;
                SingleKOResults(i).(name{counter2}). (objectives{counter}).grRateKO=grRateKO;
                SingleKOResults(i).(name{counter2}). (objectives{counter}).grRateWT=grRateWT;
            else
                disp( [objectives(counter),'is not present in model', num2str(i)] )
            end
        end
    end
end
save singleKO