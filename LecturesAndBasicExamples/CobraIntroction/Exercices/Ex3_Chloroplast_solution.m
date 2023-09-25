%% Producing the Chloroplast Carbon Metabolism Model directly in the COBRA toolbox
clear all;clc;

%change the solver
solverOK=changeCobraSolver('ibm_cplex');
%% a) Create a constraint-based model of a simple biochemical system directly in COBRA using the createModel command.
% Creat a new COBRA model with the following reactions

ReactionFormulas = {'RuBP + CO2 -> 2 PGA','PGA + ATP <=> BPGA',...
    'BPGA + NADPH <=> GAP','FBP -> F6P','SBP -> S7P',...
    'Ru5P + ATP -> RuBP','G1P + ATP ->',...
    'ATP <=>','GAP <=> DHAP','DHAP + GAP <=> FBP','GAP + F6P <=> E4P + X5P',...
    'DHAP + E4P <=> SBP','GAP + S7P <=> R5P + X5P','R5P <=> Ru5P',...
    'X5P <=> Ru5P','F6P <=> G6P','G6P <=> G1P','-> G1P','GAP ->',...
    'G6P -> R5P + 2 NADPH + CO2','NADPH <=>','CO2 <=>',...
    'F6P + E4P <=> GAP + S7P'};

ReactionNames = {'v1','v2','v3','v4','v5','v6','v7','v8','v9','v10',...
    'v11','v12','v13','v14','v15','v16','v17','v18','v19','v20','v21','v22','v23'};
GeneNames={'Gene1','Gene2','Gene3','Gene4','Gene5','Gene6','Gene7','Gene8','Gene9','Gene10',...
    'Gene11','Gene12','Gene13','Gene14','Gene15','Gene16','Gene17','Gene18','Gene19','Gene20','Gene21','Gene22','Gene23'};

orig_model = createModel(ReactionNames, ReactionNames, ReactionFormulas,'grRuleList',GeneNames);

model=orig_model;
% finally print reactions in the model
disp(' ')
disp('model contains now the following reactions:')
printRxnFormula(model,model.rxns);
disp(' ')
disp('boundaries:')
disp([model.rxns num2cell(model.lb) num2cell(model.ub)])

%% Define system boundaries: night
% RUBISCO = 0
model = changeRxnBounds(model,{'v1'},[0],'b');
% starch synthesis = 0
model = changeRxnBounds(model,{'v7'},[0],'b');
% starch degradation = 1
model = changeRxnBounds(model,{'v18'},[1],'b');
disp(' ')
disp('boundaries:')
disp([model.rxns num2cell(model.lb) num2cell(model.ub)]);

%%
% Optimization
model = changeObjective(model,'v19',1);
FBAsolution = optimizeCbModel(model);
disp(' ')
disp('FBA solution:')
printFluxVector(model,FBAsolution.x);
% [Involved_mets, Dead_ends] = draw_by_rxn(model, model.rxns, 'true', 'struc', {''}, {''}, FBAsolution.x)

disp('Night cycle')
disp('GAP production')
disp(FBAsolution.obj)%1

disp('CO2 production')
disp(FBAsolution.x(22))%3
%% Define system boundaries: day
model=orig_model;
% RUBISCO = 1
model = changeRxnBounds(model,{'v1'},[12],'b');
% starch synthesis
model = changeRxnBounds(model,{'v7'},[1],'b');
% starch degradation = 0
model = changeRxnBounds(model,{'v18'},[0],'b');
disp(' ')
disp('boundaries:')
% disp([model.rxns num2cell(model.lb) num2cell(model.ub)]);

%
% Optimization
model = changeObjective(model,'v19',1);
FBAsolution = optimizeCbModel(model);
disp(' ')
disp('FBA solution:')
printFluxVector(model,FBAsolution.x);
% [Involved_mets, Dead_ends] = draw_by_rxn(model, model.rxns, 'true', 'struc', {''}, {''}, FBAsolution.x)

disp('Day cycle')
disp('GAP production')
disp(FBAsolution.obj)%2

disp('CO2 consumption')
disp(FBAsolution.x(22))%-12

%% Gene deletion study
% single gene deletion study
[grRatio,grRateKO,grRateWT] = singleGeneDeletion(model,'FBA')
disp([model.genes num2cell(grRatio)])
%% double gene deletion study
[grRatio,grRateKO,grRateWT] = doubleGeneDeletion(model,'FBA');
disp(['_' model.genes'; model.genes num2cell(grRatio)])
