% Change the solver
% clear all
% solverOK=changeCobraSolver('ibm_cplex','all');
%% Create a new COBRA model with the following reactions

ReactionFormulas = {'-> A',' -> B','A -> C','B -> C','B -> D','C ->','D ->'};
ReactionNames = {'v1','v2','v3','v4','v5','v6','v7'};
GeneNames={'Gene1','Gene2','Gene3','Gene4','Gene5','Gene6','Gene7'};

orig_model = createModel(ReactionNames, ReactionNames, ReactionFormulas,'grRuleList',GeneNames);

disp('S'); full(orig_model.S)
disp('metabolites'); orig_model.mets
disp('reactions'); orig_model.rxns
disp(' ') 
disp('boundaries:') 
disp([orig_model.rxns num2cell(orig_model.lb) num2cell(orig_model.ub)])
disp('Genes:') 
disp([orig_model.genes])
%%
model=orig_model;
% Define system boundaries
model = changeRxnBounds(model,{'v1','v2'},[1 1],'b'); 
% finally print reactions in the model
disp(' ')
disp('model contains now the following reactions:')
printRxnFormula(model,model.rxns);
disp(' ')
disp('Boundaries:')
disp([model.rxns num2cell(model.lb) num2cell(model.ub)])
%% Optimization

model = changeObjective(model,'v6',1);
model.c
solution = optimizeCbModel(model);
disp(' ')
disp('FBA solution:')
printFluxVector(model,solution.x);

%Vizualize the flux solution
[Involved_mets, Dead_ends] = draw_by_rxn(model, model.rxns, 'true', 'struc', {''}, {''}, solution.x);

%Interactive flux analysis
surfNet(model, 'A[c]', 0, solution.x)
%% Flux Variability Analysis (FVA)

model = changeObjective(model,'v6',0)
disp('Flux Variability Analysis (FVA) with v1=[0..1]')
model = changeRxnBounds(model,{'v1'},[0],'l');%reset lower bounds to zero 

[minFlux,maxFlux] = fluxVariability(model,100);
[model.rxns num2cell(minFlux) num2cell(maxFlux)]

disp('Flux Variability Analysis (FVA) with v1=1')
model = changeRxnBounds(model,{'v1'},[1],'b');
[minFlux,maxFlux] = fluxVariability(model,100);
[model.rxns num2cell(minFlux) num2cell(maxFlux)]
%% Robustness analysis

figure
model = changeObjective(model,'v6',1); %optimize
subplot(1,2,1)
robustnessAnalysis(model,{'v4'},10)
subplot(1,2,2)
robustnessAnalysis(model,{'v7'},10)
%% Gene deletion study
% single gene deletion study

model = changeRxnBounds(model,{'v1'},[0],'l'); %reset lower bound to zero, the bounds are between 0-1
model = changeRxnBounds(model,{'v1'},[1],'u'); %reset lower bound to zero, the bounds are between 0-1

model = changeObjective(model,'v6',1); 


disp([model.rxns num2cell(model.lb) num2cell(model.ub)])

[grRatio,grRateKO,grRateWT] = singleGeneDeletion(model,'FBA')

% double gene deletion study
[grRatio,grRateKO,grRateWT] = doubleGeneDeletion(model,'FBA')
disp(['_' model.rxns'; model.rxns num2cell(grRateKO)])
%% sdsd
% dsddssd

figure
plot(1,1,'*')
% Discussion: dffsfdgsd

