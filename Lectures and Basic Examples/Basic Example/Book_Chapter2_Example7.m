%% Workbook, Chapter 2, Example 7
% Model
% uncomment the solver that you want to use
%changeCobraSolver('glpk')
changeCobraSolver('ibm_cplex')

ReactionFormulas = {'->  A','A -> B','A + B ->'};
ReactionNames = {'v1','v2','v3'};
GeneNames={'Gene1','Gene2','Gene3'};

orig_model = createModel(ReactionNames, ReactionNames, ReactionFormulas,'grRuleList',GeneNames);

%% Display and check the model structure and content e.g. via
model=orig_model;

disp('S'); full(model.S)
disp('metabolites'); model.mets
disp('reactions'); disp(model.rxns');
disp('boundaries:') 
disp([model.rxns num2cell(model.lb) num2cell(model.ub)])
disp('Genes:') 
disp([model.genes])
printRxnFormula(model)
disp('Model contains the following reactions:')
printRxnFormula(model,model.rxns);

%% Solve for a fixed v1=1
model=orig_model;

% Change reaction bounds/Define system boundaries
model = changeRxnBounds(model,{'v1'},1,'b');
disp([model.rxns num2cell(model.lb) num2cell(model.ub)])

% solve for v2 and v3
FBAsolution = optimizeCbModel(model);
disp('FBA solution optimizing:')
printFluxVector(model,FBAsolution.x);

%% Solve for a fixed v1=0.7
model=orig_model;

% Change reaction bounds/Define system boundaries
model = changeRxnBounds(model,{'v1'},0.7,'b');
disp([model.rxns num2cell(model.lb) num2cell(model.ub)])

% solve for v2 and v3
FBAsolution = optimizeCbModel(model);
disp('FBA solution optimizing:')
printFluxVector(model,FBAsolution.x);

%% Sampling for different v1 = [0 .. 1]
model=orig_model;

model = changeRxnBounds(model,{'v1'},0,'l');
model = changeRxnBounds(model,{'v1'},1,'u');
disp([model.rxns num2cell(model.lb) num2cell(model.ub)])

%Setting sampling options
options=[];
options.nWarmupPoints=100;
options.nStepsPerPoint = 1*size(model.S, 2);
options.nPointsReturned = 20;
options.nFiles=3;
[modelSampling,samples] = sampleCbModel(model,'sampleFile','ACHR',options);

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
