%% Workbook, Chapter 2, Example 15
% Model
clc
ReactionFormulas = {'->  A',' -> B','P ->','E ->','A -> B','A -> C','A -> D','B <=> C','B -> P','C + D -> E + P'};
ReactionNames = {'v1','v2','v3','v4','v5','v6','v7','v8','v9','v10'};
GeneNames={'Gene1','Gene2','Gene3','Gene4','Gene5','Gene6','Gene7','Gene8','Gene9','Gene10'};

orig_model = createModel(ReactionNames, ReactionNames, ReactionFormulas,'grRuleList',GeneNames);

%% Display and check the model structure and content e.g. via
model=orig_model;

disp('S'); full(model.S)
disp('metabolites'); model.mets
disp('reactions'); disp(model.rxns')
disp('boundaries:')
disp([model.rxns num2cell(model.lb) num2cell(model.ub)])
disp('Genes:')
disp([model.genes])

disp('Model contains the following reactions:')
printRxnFormula(model,model.rxns);

%% One solution (FBA) with blocked B uptakes while optimizing for P export
model=orig_model;

model = changeRxnBounds(model,{'v1'},1,'u');
model = changeRxnBounds(model,{'v2'},0,'b');
model = changeObjective(model,'v3',1);
disp([model.rxns num2cell(model.lb) num2cell(model.ub)])

FBAsolution = optimizeCbModel(model);
disp('FBA solution optimizing:')
printFluxVector(model,FBAsolution.x);

%% Solution overview (FVA) with blocked B uptakes while optimizing for P export
model=orig_model;

model = changeRxnBounds(model,{'v1'},1,'u');
model = changeRxnBounds(model,{'v2'},0,'b');
model = changeObjective(model,'v3',1);
disp([model.rxns num2cell(model.lb) num2cell(model.ub)])

disp('Flux Variability Analysis (FVA) with v1=1 and v4 -> max')
[minFlux,maxFlux] = fluxVariability(model,100);
disp([{'rxns' 'minFLux' 'maxFlux'}; model.rxns num2cell(minFlux) num2cell(maxFlux)]);

