%% Toy model 1 
clear all
close all

%change the solver
solverOK=changeCobraSolver('ibm_cplex','all');

%% a) Create a constraint-based model of a simple biochemical system directly in COBRA using the createModel command.
ReactionFormulas = {'->  A','A -> B','B -> C','A -> D','C -> D','D ->'};
ReactionNames = {'v1','v2','v3','v4','v5','v6'};
GeneNames={'Gene1','Gene2','Gene3','Gene4','Gene5','Gene6'};

orig_model = createModel(ReactionNames, ReactionNames, ReactionFormulas,'grRuleList',GeneNames);

%% You can display and check the model structure and content e.g. via
% Show S, metabolites, reactions, lower bounds, upper bounds, and genes
disp('S'); full(orig_model.S)
disp('metabolites'); orig_model.mets
disp('reactions'); orig_model.rxns'
disp(' ') 
disp('boundaries:') 
disp([orig_model.rxns num2cell(orig_model.lb) num2cell(orig_model.ub)])
disp('Genes:') 
disp([orig_model.genes])
%% You can print out the reaction formulas using the printRxnFormula command (check the COBRA's protocol for details).
% finally print reactions in the model
disp(' ')
disp('Model contains the following reactions:')
printRxnFormula(orig_model,orig_model.rxns);

%% b)  Fix v1 = 2 and optimize for v2, v3, or v6: What results do you get? Why?
model=orig_model;

% Change reaction bounds/Define system boundaries
model = changeRxnBounds(model,{'v1'},[2],'b');
disp([model.rxns num2cell(model.lb) num2cell(model.ub)])

%Optimization for v2:
model = changeObjective(model,'v2',1);
FBAsolution = optimizeCbModel(model);
disp(' ')
disp('FBA solution optimizing for v2:')
printFluxVector(model,FBAsolution.x);
%Vizualize the flux solution
%[Involved_mets, Dead_ends] = draw_by_rxn(model, model.rxns, 'true', 'struc', {''}, {''}, FBAsolution.x);

% Optimization for v3:
model = changeObjective(model,'v3',1);
FBAsolution = optimizeCbModel(model);
disp(' ')
disp('FBA solution optimizing for v3:')
printFluxVector(model,FBAsolution.x);
%Vizualize the flux solution
%[Involved_mets, Dead_ends] = draw_by_rxn(model, model.rxns, 'true', 'struc', {''}, {''}, FBAsolution.x);

% Optimization for v6:
model = changeObjective(model,'v6',1);
FBAsolution_sampling = optimizeCbModel(model);
disp(' ')
disp('FBA solution optimizing for v6:')
printFluxVector(model,FBAsolution_sampling.x);
%Vizualize the flux solution
%[Involved_mets, Dead_ends] = draw_by_rxn(model, model.rxns, 'true', 'struc', {''}, {''}, FBAsolution.x);
%% c)	Perform a robustness analysis for v4 as changing variable, with v3 and v5 as objective function. What is the difference? Why?

%% robustness analysis for v3 (obj) with v2 as changing variable
model = changeObjective(model,{'v3'},1);
figure; subplot(1,2,1)
%robustnessAnalysis(model, controlRxn, nPoints)
robustnessAnalysis(model, 'v2', 20)
title('Robustness analysis for v3 with v2 as changing variable')
%% robustness analysis for v4 (obj) with v2 as changing variable
model = changeObjective(model,{'v4'},1);
subplot(1,2,2)
robustnessAnalysis(model, 'v2', 10)
title('Robustness analysis for v4 with v2 as changing variable')

%% d) Perform a flux variability analysis for v1 = 2 and objective function v3, v4, or v6. 
model=orig_model;
% Define system boundaries
model = changeRxnBounds(model,{'v1'},[2],'b');

%Optimization for v3:
model = changeObjective(model,'v3',1);
% Flux Variability Analysis (FVA)
disp('Flux Variability Analysis (FVA) with v1=2 and optimizing v3')
[minFlux,maxFlux] = fluxVariability(model,100);
[{'rxns' 'minFLux' 'maxFlux'}; model.rxns num2cell(minFlux) num2cell(maxFlux)]

% Optimization for v4:
model = changeObjective(model,'v4',1);
% Flux Variability Analysis (FVA)
disp('Flux Variability Analysis (FVA) with v1=2 and optimizing v4')
[minFlux,maxFlux] = fluxVariability(model,100);
[model.rxns num2cell(minFlux) num2cell(maxFlux)]

% Optimization for v6:
model = changeObjective(model,'v6',1);
% Flux Variability Analysis (FVA)
disp('Flux Variability Analysis (FVA) with v1=2 and optimizing v6')
[minFlux,maxFlux] = fluxVariability(model,100);
[model.rxns num2cell(minFlux) num2cell(maxFlux)]

%% Perform a flux variability analysis for v1 = v4 = 2 and objective function v5. Discuss the difference.
model = changeRxnBounds(model,{'v1','v4'},[2],'b');
model = changeObjective(model,'v6',1);
disp('Flux Variability Analysis (FVA) with v1=v4=2')
[minFlux,maxFlux] = fluxVariability(model,100);
[model.rxns num2cell(minFlux) num2cell(maxFlux)]

%% e) Perform a sampling of allowed flux distributions for v1 = 2. Discuss the result.

model=orig_model;

% Define system boundaries
model = changeRxnBounds(model,{'v1'},[2],'b');

model = changeRxnBounds(model,'v6',0.9*FBAsolution_sampling.f,'l');

%Setting sampling options
options.nStepsPerPoint = 8 * size(model.S, 2);
options.nPointsReturned = 1000;

%Sampling
tic
[modelSampling,samples] = sampleCbModel(model,'sampleFile','ACHR',options);
toc

%Visualization of the sampling results:
% sampleScatterMatrix(model.rxns, model,samples, [], [], true) 
%sampleScatterMatrix({'v2','v3','v4'},modelSampling,samples)
%% Visualization of the sampling results (correlation between reactions)
figure
for counter=1:size(samples,1)
    for counter2=counter:size(samples,1)
        if counter~=counter2
            subplot(size(samples,1),size(samples,1),(counter-1)*size(samples,1)+counter2)
            plot(samples(counter,:),samples(counter2,:),'*')
            xlabel(model.rxns(counter)),ylabel(model.rxns(counter2))
        end
    end
end

%% Visualization of the distribution of the sampling results
% rxnsIdx=1:numel(model.rxns);
% figure
% for i = rxnsIdx
% nbins = 20;
% [yUn, xUn] = hist(samples(i, :), nbins);
% subplot(3, 2, find(rxnsIdx==i))
% h1 = plot(xUn, yUn);
% xlabel('Flux (mmol/gDW/h)')
% ylabel('# samples')
% title(sprintf('%s (%s)', strjoin(model.subSystems{i},';'), model.rxns{i}), 'FontWeight', 'normal')
% ylim = get(gca, 'ylim');
% end

% or simply:
plotSampleHist(modelSampling.rxns, {samples}, modelSampling,[],[2,3]);

%% Identify sets of correlated reactions from the sampling data.
[sets,setNumber,setSize]=identifyCorrelSets(model,samples,1-1e-8);
for counter=1:size(sets)
    temp=cell2mat(sets(counter));
    disp(['Set ' num2str(counter) ':'])
    disp([num2cell(temp.set) temp.names])
end

%% f) Perform a single gene deletion study with the objective function v6. Discuss the result.
model=orig_model;
model = changeObjective(model,{'v6'},1);
% Single gene deletion study
[grRatio,grRateKO,grRateWT] = singleGeneDeletion(model,'FBA')

%% g) Perform a double gene deletion study. What is the only non-lethal double deletion?  
% Double gene deletion study
[grRatio,grRateKO,grRateWT] = doubleGeneDeletion(model,'FBA');
disp(['_' model.genes'; model.genes num2cell(grRatio)])
