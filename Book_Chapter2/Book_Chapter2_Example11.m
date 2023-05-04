%% Workbook, Chapter 2, Example 11
% Model
clc
ReactionFormulas = {'->  G6P','G6P ->','G6P -> 2 DHAP','DHAP -> PEP','PEP ->'};
ReactionNames = {'v1','v2','v3','v4','v5'};
GeneNames={'Gene1','Gene2','Gene3','Gene4','Gene5'};

orig_model = createModel(ReactionNames, ReactionNames, ReactionFormulas,'grRuleList',GeneNames);

%% Display and check the model structure and content e.g. via
model=orig_model;

disp('S'); full(model.S)
disp('metabolites'); model.mets
disp('reactions'); model.rxns'
disp('boundaries:')
disp([model.rxns num2cell(model.lb) num2cell(model.ub)])
disp('Genes:')
disp([model.genes])

disp('Model contains the following reactions:')
printRxnFormula(model,model.rxns);

%% Solve for fixed v1=1 and v5=1.4
model=orig_model;

% Change reaction bounds/Define system boundaries
model = changeRxnBounds(model,{'v1'},[1],'b');
model = changeRxnBounds(model,{'v5'},[1.4],'b');
disp([model.rxns num2cell(model.lb) num2cell(model.ub)])

% solve for v2 and v3
FBAsolution = optimizeCbModel(model);
disp('FBA solution optimizing:')
printFluxVector(model,FBAsolution.x);

%% FVA  for fixed v1=1 and v5=1.4
model=orig_model;
model = changeRxnBounds(model,{'v1'},[1],'b');
model = changeRxnBounds(model,{'v5'},[1.4],'b');

disp('Flux Variability Analysis (FVA) with v1=1 and v5=1.4')
[minFlux,maxFlux] = fluxVariability(model,100);
[{'rxns' 'minFLux' 'maxFlux'}; model.rxns num2cell(minFlux) num2cell(maxFlux)]

%% Additional: FVA  for fixed v1=1 only
model=orig_model;
model = changeRxnBounds(model,{'v1'},[1],'b');

disp('Flux Variability Analysis (FVA) with v1=1 only')
[minFlux,maxFlux] = fluxVariability(model,100);
[{'rxns' 'minFLux' 'maxFlux'}; model.rxns num2cell(minFlux) num2cell(maxFlux)]
