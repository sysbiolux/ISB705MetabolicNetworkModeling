clear all

%In this script we want to enhable the model to produce Luteine and S_epsilon_45_Carotene
% two demand reaction were added : epsilon carotene harvesting Luteinharvesting

%set the path
path='C:\Users\maria.pacheco\Documents\GitHub\basicAnalysis';
load(strcat(path,'\Model\Arabidopsis.mat'), 'model2010');
changeCobraSolver('ibm_cplex');
% check if the model can produce Lutein
model2010=changeObjective(model2010,'Luteinharvesting');
sol_ori=optimizeCbModel(model2010);
sol_ori.f

% add demands to deadend
[ model_corrected,deadEnd]=findDeadEndsFastbox(model2010)
model_corrected.rxnNames(end+1:numel(model_corrected.rxns))=model_corrected.rxns(numel(model_corrected.rxnNames)+1:numel(model_corrected.rxns));
model_corrected.c(end+1:numel(model_corrected.rxns))=0;
model=model_corrected;
model.c=zeros(numel(model.rxns),1);
% try to optimize again
model=changeObjective(model,'Luteinharvesting');
C=find(contains(model.rxns, 'Luteinharvesting'));
artificial=setdiff(model_corrected.rxns, model2010.rxns),
t=find(ismember(model2010.rxns,model_corrected.rxns)); %indices of medium exchange reactions
numel(t)
model.S(:,C)=model.S(:,C)*1000;
A=fastcore_4_rfastcormics(C, model, 1e-4, t);

metabolites={'Lycopene', 'Phytophoeme', 'Geranyl', 'S__40_R_41__45_Mevalonate[c]', 'S_Pyruvate[c]', 'S_Glucose[c]'};
mets=model.mets(contains(model.mets, metabolites));
Rxns2=findRxnsFromMets(model,mets);
Table=table(Rxns2, printRxnFormula(model, Rxns2));
Rxns={'R02245_c','R07916_c', 'RN03_c', 'R00200_c'}
to_close_keep=[];
to_close=[];

% check if the reactions are unblocked
A_tmp=fastcc_4_rfastcormics(model,1e-4,1)
model_temp=removeRxns(model, model.rxns(setdiff(1:numel(model.rxns),A_tmp)))
setdiff(Rxns, model.rxns(A_tmp))
% remove demand rxns that are not required for the reactions of interest to
% carry a flux
for counter =1: numel(Rxns)    
    if isempty(intersect(Rxns(counter), model.rxns(A)))
        tmp=model;
        tmp.lb(ismember(model.rxns,Rxns(counter)))=1000
        sol=optimizeCbModel(tmp)
        sol.f
        if sol.f~=0 && ~isempty(sol.f) && ~isnan(sol.f)
            A1=fastcore_4_rfastcormics(C, tmp, 1e-4, t);
            to_close=setdiff(artificial, model.rxns(A1));
            to_close_keep=[to_close_keep;to_close];
            'close'
            numel(unique(to_close_keep))
            model.lb(ismember(model.rxns, to_close_keep))=0;
            model.ub(ismember(model.rxns, to_close_keep))=0;
            A2=fastcore_4_rfastcormics(C, model, 1e-4, t);
            to_close=setdiff(artificial, model.rxns(A2));
            to_close_keep=unique([to_close_keep;to_close]);
            A=A2

        else
            
        end
        
        while  isempty(intersect(model.rxns(A2), Rxns(counter))) & numel(unique(to_close))< numel(unique(to_close_keep))
            if isempty(intersect(model.rxns(A2), Rxns(counter)))
                tmp=model;
                tmp.lb(ismember(model.rxns,Rxns(counter)))=100;
                sol=optimizeCbModel(tmp);
                sol.f
                if sol.f~=0
                    A=fastcore_4_rfastcormics(C, tmp, 1e-4, t);
                    to_close=setdiff(artificial, model.rxns(A));
                    to_close_keep=unique([to_close_keep;to_close]);
                  
                    model.lb(ismember(model.rxns, to_close_keep))=0;
                    model.ub(ismember(model.rxns, to_close_keep))=0;
                    A=fastcore_4_rfastcormics(C, model, 1e-4, t);
                    
                else
                    sss
                end
                
                
            else
                vvv
            end
            
        end
        
    end
end

Table=table(model.rxns(A),printRxnFormula(model,model.rxns(A)))
sol=optimizeCbModel(model, 'max', 'zero')
table(model.rxns(sol.x~=0),printRxnFormula(model,model.rxns(sol.x~=0)))

model=addReaction(model, 'diphosphomevalonate decarboxylase','S__40_R_41__45_5_45_Diphosphomevalonate[c]  + S_ATP[c] -> S_ADP[c] + S_Phosphatidate[c] + S_Isopentenyl_32_diphosphate[c]   + S_CO2[c]') 
model=addReaction(model, 'CO2_trans','S_CO2[c] -> S_CO2[e]') 
model=addReaction(model, 'CO2_ex','S_CO2[e] -> ') 
model=addReaction(model, 'PI_trans','S_Phosphatidate[c] -> S_Phosphatidate[e]') 
model=addReaction(model, 'S_Phosphatidate_ex','S_Phosphatidate[e] -> ') 
model=removeRxns(model, {'SinkToRemoveDeadEnd_I_S_1_45__40_2_45_Carboxyphenylamino_41__45_1_39__45_deoxy_45_D_45_ribulose_32_5_39__45_phosphate[c]_S_1_45__40_2_45_Carboxyphenylamino_41__45_1_39__45_deoxy_45_D_45_ribulose_32_5_39__45_phosphate[c]'})

Rxns=union(Rxns,model.rxns(end-4:end))
%%
for counter =1: numel(Rxns)
    if isempty(intersect(Rxns(counter), model.rxns(A)))
        tmp=model;
        tmp.lb(ismember(model.rxns,Rxns(counter)))=1000;
        sol=optimizeCbModel(tmp)
        sol.f
        if sol.f~=0 && ~isempty(sol.f) & ~isnan(sol.f)
            A=fastcore_4_rfastcormics(C, tmp, 1e-4, t);
            to_close=setdiff(artificial, model.rxns(A));
            to_close_keep=[to_close_keep;to_close];
            'close'
            numel(unique(to_close_keep))
            model.lb(ismember(model.rxns, to_close_keep))=0;
            model.ub(ismember(model.rxns, to_close_keep))=0;
            A=fastcore_4_rfastcormics(C, model, 1e-4, t);
            
        else
        end
        while isempty(intersect(model.rxns(A), Rxns(counter))) & numel(unique(to_close))< numel(unique(to_close_keep))
         
            
            
            if isempty(intersect(model.rxns(A), Rxns(counter)))
                tmp=model;
                tmp.lb(ismember(model.rxns,Rxns(counter)))=100;
                sol=optimizeCbModel(tmp);
                sol.f
                if sol.f~=0
                    A=fastcore_4_rfastcormics(C, tmp, 1e-4, t);
                    to_close=setdiff(artificial, model.rxns(A));
                    to_close_keep=[to_close_keep;to_close];
                    'close'
                    numel(unique(to_close_keep))
                    model.lb(ismember(model.rxns, to_close_keep))=0;
                    model.ub(ismember(model.rxns, to_close_keep))=0;
                    A=fastcore_4_rfastcormics(C, model, 1e-4, t);
                    
                else
                    sss
                end
                
                
            else
                vvv
            end
            
        end
        
    end
end
Table=table(model.rxns(A),printRxnFormula(model,model.rxns(A)))
model=removeRxns(model,'SinkToRemoveDeadEnd_I_S_2_44_2_45_Dichloroacetaldehyde[c]_S_2_44_2_45_Dichloroacetaldehyde[c]')
model=removeRxns(model,'SinkToRemoveDeadEnd_I_S_2_45_Methyl_45_1_45_hydroxybutyl_45_ThPP[c]_S_2_45_Methyl_45_1_45_hydroxybutyl_45_ThPP[c]')
model=removeRxns(model,'SinkToRemoveDeadEnd_I_S_2_44_5_45_Dihydroxypyridine[c]_S_2_44_5_45_Dihydroxypyridine[c]')

sol=optimizeCbModel(model, 'max', 'zero')
table(model.rxns(sol.x~=0),printRxnFormula(model,model.rxns(sol.x~=0)), sol.x(sol.x~=0))
