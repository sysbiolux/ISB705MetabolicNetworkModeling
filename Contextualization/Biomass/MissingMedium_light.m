% Load data for the example 
load('brain_model.mat');
load('medium_example.mat');

%% Specific fix for this example 
sum(ismember(brain_model.rxns,{'DM_4abut[c]'}))
to_remove = find(ismember(brain_model.rxns, {'DM_4abut[c]'}));
model_orig = removeRxns(brain_model, brain_model.rxns(to_remove(1))); % remove duplicated rxns

medium_exchanges = cellstr(intersect(medium_example, brain_model.mets));
epsilon = 1e-4;
model_orig.metFormulas(cellfun('isempty', model_orig.metFormulas)) = cellstr('C');
model_orig.metFormulas(strcmp(model_orig.metFormulas, 'X')) = cellstr('C');

%% UNCONSTRAINED

model = model_orig;
model = changeObjective(model, 'biomass_maintenance');
printObjective(model);
solution_orig = optimizeCbModel(model, 'max');
disp(solution_orig.f) %288

%% MEDIUM CONSTRAINED

exRxns = model.rxns(findExcRxns(model)); % All exchanges reactions
medium_exchanges = strrep(medium_exchanges, '[e]', '[csf]'); % (Specific to a model type ?)
all_medium_rxns = findRxnsFromMets(model, medium_exchanges); % Rxns of the mets in the medium 

exToKeep = intersect(exRxns, all_medium_rxns); % Medium exchange reactions (MER)
disp(numel(exToKeep)) % Number of MER

% Remove uptake of carbon sources which are not in the medium from the model
[model_cons] = constrain_model_rFASTCORMICS(model, medium_exchanges, [], 'biomass_maintenance', 'biomass_maintenance'); 

exToRemove = model.rxns(model.lb ~= model_cons.lb | model.ub ~= model_cons.ub); % Ex rxns that are close when we do the medium constrained

solution_cons = optimizeCbModel(model_cons, 'max');
disp(strcat('Biomass production after medium constraints: ', num2str(solution_cons.f)));

%% CONSISTENT MODEL 

model = model_orig;
model = changeObjective(model, 'biomass_maintenance');

% Find missing metabolites 
C = find(ismember(model.rxns, 'biomass_maintenance')); % Index of the objective rxn 

allowedInputRxns = find(ismember(model.rxns, exToKeep)); % Indices of MER
disp(strcat('Allowed input rxns: ', num2str(numel(allowedInputRxns)))); % Number of MER 

model.S(:,C)=model.S(:,C)*1000; %TS: multiple biomass coefficients to overcome numerical issues

% Using fastcore_4_rfastcormics 
% INPUT :     
    % C = core rxns -> objective rxn 
    % f = unpenalized rxns 
% OUTPUT : 
    % A = List of rxns index that should be kept 
 
A = fastcore_4_rfastcormics(C, model, epsilon, allowedInputRxns); 
% rxns_to_keep = model.rxns(A); % List of the minimal rxns to keep the objective function 

disp(strcat('Reactions in minimal mode to produce biomass:', num2str(numel(A))));

% Check if the objective function is in ('biomass_maintenance')
if ~isempty(setdiff(C,A))
    error('Objective function is missing')
else 
    disp('Objective function is present')
end

% Missing rxns in the medium : 
neededUptakes = intersect(exToRemove, model.rxns(A));

%% APPROACH 1 : MEDIUM CONSTRAINED CONSISTENT MODEL 
% Time : 0.222 s

model = model_orig;
model = changeObjective(model, 'biomass_maintenance');

% STEP 1 : Consistent model (+ allowed input)
model_min = removeRxns(model, model.rxns(setdiff(1:numel(model.rxns), A)));

% STEP 2 : Add medium constrained to the consistent model
% Remove uptake of carbon sources which are not in the medium from the model
[model_cons] = constrain_model_rFASTCORMICS(model, medium_exchanges, [], 'biomass_maintenance', 'biomass_maintenance');
% Reset the bounds of the needed uptakes rxns
model_cons.lb(ismember(model.rxns, neededUptakes)) = model_orig.lb(ismember(model.rxns, neededUptakes)); 
model_cons.ub(ismember(model.rxns, neededUptakes)) = model_orig.ub(ismember(model.rxns, neededUptakes));

solution = optimizeCbModel(model_cons,'max');
solution.f

%% APPROACH 2 : MEDIUM CONSTRAINED CONSISTENT MODEL
% Time : 28.341 s (fastcc_4_rfastcormics)

model = model_orig;
model = changeObjective(model, 'biomass_maintenance');

% add needed metabolites to the medium
missing_needed_metabolites = findMetsFromRxns(model, neededUptakes)
medium_with_missing_metabolites = cat(1, medium_exchanges, missing_needed_metabolites);

[EX] = findExcRxns(model); % Find all exchange reactions
[Ex_open] = findRxnsFromMets(model, medium_with_missing_metabolites);
Ex_to_close = setdiff(model.rxns(EX),Ex_open);

medium_model = model; % Save in case of mistake
medium_model.lb(ismember(medium_model.rxns, Ex_to_close)) = 0

A = fastcc_4_rfastcormics(medium_model, 1e-4, 0);
medium_constrained_consistent_model = removeRxns(medium_model, medium_model.rxns(setdiff(1:numel(medium_model.rxns),A)))

solution = optimizeCbModel(medium_constrained_consistent_model,'max');
solution.f
