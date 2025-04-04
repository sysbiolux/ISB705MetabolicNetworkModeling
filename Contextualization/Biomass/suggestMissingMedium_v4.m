%% load model & medium
clearvars -except solverOK, clc, close all
% addpath('C:\Users\maria.pacheco\OneDrive - University of Luxembourg\Documents\GitHub\rFASTCORMICS')
% addpath('C:\Users\maria.pacheco\Documents\GitHub\rFASTCORMICS')
addpath(genpath('C:\Users\thomas.sauter\OneDrive - University of Luxembourg\uni.lu\MISB\ISB705\SelfStudyCourse_2023\ISB705MetabolicNetworkModeling-main'))
% changeCobraSolver('ibm_cplex');

load('brain_model.mat');
load('medium_example.mat');

%% fix a reaction which is present twice in the model
sum(ismember(brain_model.rxns,{'DM_4abut[c]'}))
to_remove=find(ismember(brain_model.rxns,{'DM_4abut[c]'}));
model_orig=removeRxns(brain_model, brain_model.rxns(to_remove(1))); % remove duplicated rxns

medium_exchanges = cellstr(intersect(medium_example, brain_model.mets));
epsilon=1e-4;

%% find metabolites without formulas (and add C there)
model_orig.metFormulas(cellfun('isempty',model_orig.metFormulas))=cellstr('C');
model_orig.metFormulas(strcmp(model_orig.metFormulas,'X'))=cellstr('C');

%% biomass of unconstrained model
model=model_orig;
model = changeObjective(model,'biomass_maintenance');
checkObjective(model);
solution_orig = optimizeCbModel(model,'max');
disp(solution_orig.f) %288

%% biomass of medium comstrained model
exRxns=model.rxns(findExcRxns(model));
medium_exchanges=strrep(medium_exchanges,'[e]','[csf]');
temp=findRxnsFromMets(model,medium_exchanges);
% convert
exToKeep=intersect(exRxns,temp); %medium exchange reactions
disp(numel(exToKeep))

% remove uptake of carbon sources which are not in the medium from the model
[model_cons] = constrain_model_rFASTCORMICS(model,medium_exchanges , [], 'biomass_maintenance', 'biomass_maintenance');
exToRemove=model.rxns(model.lb~=model_cons.lb | model.ub~= model_cons.ub);

solution_cons = optimizeCbModel(model_cons,'max');
disp(strcat(' biomass production after medium constraints: ',num2str(solution_cons.f)));

%% determine missing metabolites uptakes
model=model_orig;
model = changeObjective(model,'biomass_maintenance');
C=find(ismember(model.rxns,'biomass_maintenance'));
allowedInputRxns=find(ismember(model.rxns,exToKeep)); %indices of medium exchange reactions
disp(' ')
disp(strcat('allowed input rxns: ', num2str(numel(allowedInputRxns))));
model.S(:,C)=model.S(:,C)*1000; %TS: multiple biomass coefficients to overcome numerical issues

% run fastcore for biomass only, with medium exchange reactions unpenalized
A=fastcore_4_rfastcormics(C, model, epsilon, allowedInputRxns);
model.rxns(A);
disp(' ')
disp('Reactions in minimal mode to produce biomass:')
numel(A)
% check if Biomass is in
if ~isempty(setdiff(C,A))
    error('No biomass produced')
end
disp('Missing medium metabolites (substrates & products):')
neededUptakes=intersect(exToRemove, model.rxns(A))

%% check biomass of minimal model with A only (including additional metabolites)
model=model_orig;

model = changeObjective(model,'biomass_maintenance');
model_min = removeRxns(model,model.rxns(setdiff(1:numel(model.rxns),A))); %keep only A
% Acc =  fastcc_4_rfastcormics(model_min,epsilon,0);
% numel(Acc)
[~,~,IB]=intersect(model_min.rxns,exRxns);
exRxn_minimal=exRxns(IB);
solution_min = optimizeCbModel(model_min,'max');
disp(strcat('minimal solution: ', num2str(solution_min.f)))

[~, IA, IB]=intersect(model_min.rxns,exRxn_minimal); % all exchanges
T=table(exRxn_minimal(IB),solution_min.x(IA));

% defined as export in S matrix
[~,exports]=find(model_orig.S(:,ismember(model_orig.rxns, neededUptakes))<0);
exports=neededUptakes(exports);

% defined as import in S matrix
[~,imports]=find(model_orig.S(:,ismember(model_orig.rxns, neededUptakes))>0);
imports=neededUptakes(imports);

% keep only importing new medium components (depending in export/import
% definition)
[~, IA, IB]=intersect(model_min.rxns, exports);
exports=exports(IB(solution_min.x(IA)<0));
[I, IA, IB]=intersect(model_min.rxns, exports);
imports=imports(IB(solution_min.x(IA)>0));

disp('Missing medium metabolites (substrates only):')
neededUptakes=union(imports, exports)

%% check biomass of minimal model with medium and additional metabolites
model=model_orig;
model = changeObjective(model,'biomass_maintenance');
[model_cons] = constrain_model_rFASTCORMICS(model,medium_exchanges , [], 'biomass_maintenance', 'biomass_maintenance');
model_cons.lb(ismember(model.rxns,exToKeep))=model_orig.lb(ismember(model.rxns,exToKeep));
model_cons.ub(ismember(model.rxns,exToKeep))=model_orig.ub(ismember(model.rxns,exToKeep));

model_cons.lb(ismember(model.rxns,neededUptakes))=model_orig.lb(ismember(model.rxns,neededUptakes));
model_cons.ub(ismember(model.rxns,neededUptakes))=model_orig.ub(ismember(model.rxns,neededUptakes));

solution_corr = optimizeCbModel(model_cons,'max');
solution_corr.f

% determine minimal uptakes of additiional medium components which still
% gives biomass and uses Glc (quick & dirty)
if abs(solution_corr.x((ismember(model_cons.rxns, 'EX_glc_D[csf]'))))<epsilon && solution_corr.f>epsilon
    solution_corr.x((ismember(model_cons.rxns, neededUptakes)))
    sol=solution_corr;
    tmp=model_cons;
    while sol.f>0 % && abs(sol.x(ismember(model_cons.rxns, 'EX_glc_D[csf]')))<epsilon
        tmp.lb(ismember(model_cons.rxns, neededUptakes))=tmp.lb(ismember(model_cons.rxns, neededUptakes))+10; %step-size
        sol = optimizeCbModel(tmp,'max');
        sol.f
        if  sol.f >epsilon  && abs(sol.x(ismember(model_cons.rxns, 'EX_glc_D[csf]')))>epsilon
            final=tmp;
        end
    end
    sol = optimizeCbModel(final,'max'); %118
end

[~, IA, IB]=intersect(final.rxns, exRxns);
temp=table(exRxns(IB),printRxnFormula(final,exRxns(IB)),sol.x(IA));
res=temp(abs(temp.Var3)>epsilon,:);
res = sortrows(res,'Var3','ascend')
