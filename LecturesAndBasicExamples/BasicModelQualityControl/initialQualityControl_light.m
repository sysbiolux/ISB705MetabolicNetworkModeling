% initial quality control (light version), T.Sauter from 07/23 on
clearvars -except solverOK
clc

model=readCbModel('..\..\Model\Blautia_hydrogenotrophica_DSM_10507.xml')

%% annotation & stats
% Are the following fields existing and correct?
model.mets(1:10)
model.rxns(1:10)
model.genes(1:10)
model.geneNames(1:10)
model.rules(1:10)
model.comps
model.subSystems{1:3}

% Number of metabolites, reactions and genes
numel(model.mets)
numel(model.rxns)
numel(model.genes)

% Number of subsystems & Number of reactions per subsystem
temp=model.subSystems;
% subSystems=reshape([temp{:}],size(temp,1),size(temp{1},2));
% subSystems=[temp{:}];
subSystems=[];
for counter=1:numel(model.subSystems)
   if isempty(model.subSystems{counter})
       subSystems=[subSystems, 'NA'];
   else
       subSystems=[subSystems, temp{counter}];
   end
end
numel(unique(subSystems))
[uc, ~, idc] = unique( subSystems ) ;
counts = accumarray( idc, ones(size(idc)) ) ;
table(uc',counts)

% Values of lower and upper bounds
unique(model.lb)
unique(model.ub)

% Reactions in objective function
unique(model.c)
model.rxns(find(model.c==1))

%% flux consistency
[A, modelFlipped, V] = fastcc(model, 1e-4, 2);
% Number and ratio of consistent reactions
numel(A)
numel(A)/numel(model.rxns)

% Number of blocked/closed reactions
blocked=setdiff(1:numel(model.rxns),A);
numel(blocked)

%% Maximal growth rate
sol=optimizeCbModel(model,'max','zero') %minimal solution

% Minimal mode for enabling growth and respective number of reactions per subsystem
sum(sol.v~=0) %number of non-zero fluxes
temp=subSystems(find(sol.v~=0));
[uc, ~, idc] = unique( temp ) ;
counts = accumarray( idc, ones(size(idc)) ) ;
table(uc',counts) %involved subSystems

% Minimal medium for enabling growth
intersect(model.rxns(find(sol.v<0)),model.rxns(findExcRxns(model)))
