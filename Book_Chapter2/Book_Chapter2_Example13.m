%% Workbook, Chapter 2, Example 13
% Model
clc
ReactionFormulas = {'->  A','A ->','A ->'};
ReactionNames = {'v1','v2','v3'};
GeneNames={'Gene1','Gene2','Gene3'};

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

%% Solution space with enzyme capacity only
model=orig_model;

% Enzyme capacity v1<=1
model = changeRxnBounds(model,{'v1'},[1],'u');
disp([model.rxns num2cell(model.lb) num2cell(model.ub)])

disp('Flux Variability Analysis (FVA) with v1=1 and v5=1.4')
[minFlux,maxFlux] = fluxVariability(model,100);
[{'rxns' 'minFLux' 'maxFlux'}; model.rxns num2cell(minFlux) num2cell(maxFlux)]

%% Solution space with enzyme capacity and optimization of v2
model=orig_model;

% Enzyme capacity v1<=1
model = changeRxnBounds(model,{'v1'},[1],'u');
disp([model.rxns num2cell(model.lb) num2cell(model.ub)])

%Optimization for v2:
model = changeObjective(model,'v2',1);

disp('Flux Variability Analysis (FVA) with v1=1 and v5=1.4')
[minFlux,maxFlux] = fluxVariability(model,100);
[{'rxns' 'minFLux' 'maxFlux'}; model.rxns num2cell(minFlux) num2cell(maxFlux)]

%% Additionally v2<=0.6
model=orig_model;

% Enzyme capacity v1<=1 and v2 <=0.6
model = changeRxnBounds(model,{'v1'},[1],'u');
model = changeRxnBounds(model,{'v2'},[0.6],'u');
disp([model.rxns num2cell(model.lb) num2cell(model.ub)])

%Optimization for v2:
model = changeObjective(model,'v2',1);

disp('Flux Variability Analysis (FVA) with v1=1 and v5=1.4')
[minFlux,maxFlux] = fluxVariability(model,100);
[{'rxns' 'minFLux' 'maxFlux'}; model.rxns num2cell(minFlux) num2cell(maxFlux)]

%% What would happen if we impose v2=0.6 and v3=0.5?
model=orig_model;

% Enzyme capacity v1<=1 and v2 <=0.6
model = changeRxnBounds(model,{'v1'},[1],'u');
model = changeRxnBounds(model,{'v2'},[0.6],'b');
model = changeRxnBounds(model,{'v3'},[0.5],'b');
disp([model.rxns num2cell(model.lb) num2cell(model.ub)])

%Optimization for v2:
model = changeObjective(model,'v2',1);

disp('Flux Variability Analysis (FVA) with v1=1 and v5=1.4')
[minFlux,maxFlux] = fluxVariability(model,100);
[{'rxns' 'minFLux' 'maxFlux'}; model.rxns num2cell(minFlux) num2cell(maxFlux)]
