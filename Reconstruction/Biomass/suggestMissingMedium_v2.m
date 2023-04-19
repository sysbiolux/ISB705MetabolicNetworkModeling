%% load model & medium
% addpath('C:\Users\thomas.sauter\Documents\MATLAB\rFASTCORMICS_0320')
load('brain_model.mat')
model_orig=brain_model
epsilon = 1e-4

w = readtable("Thiele2020_CSF_.csv");
CSF_medium=table2cell(w(:,1));
z = strcat(CSF_medium, "[csf]");
u = cellstr(intersect(z, brain_model.mets))

%% biomass of unconstrained model
model=model_orig
model = changeObjective(model,'biomass_maintenance');
checkObjective(model);
solution = optimizeCbModel(model,'max');
solution.f

%% biomass of medium comstrained model
model=model_orig
model = changeObjective(model,'biomass_maintenance');

exRxns=model.rxns(findExcRxns(model));
temp=findRxnsFromMets(model,u);
exToKeep=intersect(exRxns,temp); %medium exchange reactions
numel(exToKeep)
exToRemove=setdiff(exRxns,exToKeep); %no medium exchange reactions
numel(exToRemove)

model.lb(find(ismember(model.rxns,exRxns)))=0;
model.ub(find(ismember(model.rxns,exRxns)))=0;

%keep only medium exchange reactions
model.lb(find(ismember(model.rxns,exToKeep)))=model_orig.lb(find(ismember(model.rxns,exToKeep)));
model.ub(find(ismember(model.rxns,exToKeep)))=model_orig.ub(find(ismember(model.rxns,exToKeep)));

solution = optimizeCbModel(model,'max');
solution.f

%% determine missing metabolites uptakes
model=model_orig;
model = changeObjective(model,'biomass_maintenance');
C=find(ismember(model.rxns,'biomass_maintenance'))
t=find(ismember(model.rxns,exToKeep)); %indices of medium exchange reactions
numel(t)

model.S(:,C)=model.S(:,C)*1000; %TS: multiple biomass coefficients to overcome numerical issues

% Optional: Exclude some metabolites completely, e.g. 'EX_dhap[csf]'
% model.lb(find(ismember(model.rxns,'EX_dhap[csf]')))=0;
% model.ub(find(ismember(model.rxns,'EX_dhap[csf]')))=0;

% run fastcore for biomass only, with medium exchange reactions unpenalized
A=fastcore_4_rfastcormics(C, model, epsilon, t); 
model.rxns(A);
numel(A)

% minimal set of additional non-medium metabolite (uptakes)
neededUptakes=intersect(exToRemove, model.rxns(A))
% neededUptakes=intersect(exRxns, model.rxns(A))

%% check biomass of minimal model with A only (including additional metabolites)
model=model_orig;
model = changeObjective(model,'biomass_maintenance');
model_min = removeRxns(model,model.rxns(setdiff(1:numel(model.rxns),A)));
model_min = removeRxns(model_min,setdiff(model_min.rxns,model.rxns(A)));

Acc =  fastcc_4_rfastcormics(model_min,epsilon,0);
numel(Acc)

solution = optimizeCbModel(model_min,'max');
solution.f

%% check biomass of minimal model with medium and additional metabolites
model=model_orig
model = changeObjective(model,'biomass_maintenance');

model.lb(find(ismember(model.rxns,exRxns)))=0;
model.ub(find(ismember(model.rxns,exRxns)))=0;

model.lb(find(ismember(model.rxns,exToKeep)))=model_orig.lb(find(ismember(model.rxns,exToKeep)));
model.ub(find(ismember(model.rxns,exToKeep)))=model_orig.ub(find(ismember(model.rxns,exToKeep)));

model.lb(find(ismember(model.rxns,neededUptakes)))=model_orig.lb(find(ismember(model.rxns,neededUptakes)));
model.ub(find(ismember(model.rxns,neededUptakes)))=model_orig.ub(find(ismember(model.rxns,neededUptakes)));

% model.lb(find(ismember(model.rxns,exRxns)))=model_orig.lb(find(ismember(model.rxns,exRxns)));
% model.ub(find(ismember(model.rxns,exRxns)))=model_orig.ub(find(ismember(model.rxns,exRxns)));

solution = optimizeCbModel(model,'max');
solution.f
