%% Toy model 2
clearvars -except solverOK
close all force, clc
% a) Create a constraint-based model of a simple biochemical system directly in COBRA using the createModel command.
ReactionFormulas = {'0.5 BC1 + 0.5 BC2 ->','2 A + C -> BC1','C + 3 D -> BC2',...
    ' <=> A','A -> B','B -> D','A <=> C','C -> D','D <=> C','D <=>'};
ReactionNames = {'mu','RBC1','RBC2','R1','R2','R3','R4','R5','R6','R7'};
GeneNames={'GeneMU','GeneBC1','GeneBC2','Gene0','Gene1','Gene2','Gene3','Gene4','Gene5','Gene6'};

orig_model = createModel(ReactionNames, ReactionNames, ReactionFormulas,'grRuleList',GeneNames)

%% set boundary condiations, e.g. measured fluxes
model=orig_model;
model = changeRxnBounds(model, {'R1'}, [1], 'b');
model = changeRxnBounds(model, {'R7'}, [0], 'l');

% b) FBA (Flux Balance Analysis) simulating maximal growth
model=changeObjective(model,'mu');
FBAsolution=optimizeCbModel(model);
printFluxVector(model, FBAsolution.x, false, false)
[Involved_mets, Dead_ends] = draw_by_rxn(model, model.rxns, 'true', 'struc', {''}, {''}, FBAsolution.x)

%% c) Flux variability
model=orig_model;
model = changeRxnBounds(model, {'R1'}, [1], 'b');
model = changeRxnBounds(model, {'R7'}, [0], 'l');

[minFlux1,maxFlux1] = fluxVariability(model,100);
rxnNames = model.rxns; %{'PGI','PFK','FBP','FBA','TPI','GAPD','PGK','PGM','ENO','PYK','PPS','G6PDH2r','PGL','GND','RPI','RPE','TKT1','TKT2','TALA'};
rxnID = findRxnIDs(model,rxnNames);
disp('Flux variability:')
printLabeledData(model.rxns(rxnID),[minFlux1(rxnID) maxFlux1(rxnID) maxFlux1(rxnID)-minFlux1(rxnID)],false)

%% Flux variability with R3=0
model=orig_model;
model = changeRxnBounds(model, {'R1'}, [1], 'b');
model = changeRxnBounds(model, {'R7'}, [0], 'l');

model = changeRxnBounds(model, {'R3','R6'}, [0; -1000], 'l');
model = changeRxnBounds(model, {'R3','R6'}, [0; 1000], 'u');
[minFlux2,maxFlux2] = fluxVariability(model,100);
rxnNames = model.rxns; %{'PGI','PFK','FBP','FBA','TPI','GAPD','PGK','PGM','ENO','PYK','PPS','G6PDH2r','PGL','GND','RPI','RPE','TKT1','TKT2','TALA'};
rxnID = findRxnIDs(model,rxnNames);
disp(' ')
disp('Flux variability with R3=0:')
printLabeledData(model.rxns(rxnID),[minFlux2(rxnID) maxFlux2(rxnID) maxFlux2(rxnID)-minFlux2(rxnID)],false)

% Flux variability with R3,R6=0
% % % model = changeRxnBounds(model, {'R3','R6'}, [0; 0], 'l');
% % % model = changeRxnBounds(model, {'R3','R6'}, [0; 0], 'u');
% % % [minFlux,maxFlux] = fluxVariability(model,100);
% % % rxnNames = model.rxns; %{'PGI','PFK','FBP','FBA','TPI','GAPD','PGK','PGM','ENO','PYK','PPS','G6PDH2r','PGL','GND','RPI','RPE','TKT1','TKT2','TALA'};
% % % rxnID = findRxnIDs(model,rxnNames);
% % % printLabeledData(model.rxns(rxnID),[minFlux(rxnID) maxFlux(rxnID) maxFlux(rxnID)-minFlux(rxnID)],false)

%% FVA similarity
[overallSim, rxnSim] = FVAsimilarity([minFlux1,maxFlux1], [minFlux2,maxFlux2])

printLabeledData(model.rxns(rxnID),[minFlux1,maxFlux1, minFlux2,maxFlux2, rxnSim],false)
