%% Toy model 2

%change the solver
changeCobraSolver('ibm_cplex');
%% 
% 
%% a) Create a constraint-based model of a simple biochemical system directly in COBRA using the createModel command.
% 

ReactionFormulas = {'0.5 BC1 + 0.5 BC2 ->','2 A + C -> BC1','C + 3 D -> BC2',...
    ' <=> A','A -> B','B -> D','A <=> C','C -> D','D <=> C','D <=>'};
ReactionNames = {'mu','RBC1','RBC2','R1','R2','R3','R4','R5','R6','R7'};
GeneNames={'GeneMU','GeneBC1','GeneBC2','Gene0','Gene1','Gene2','Gene3','Gene4','Gene5','Gene6'};

orig_model = createModel(ReactionNames, ReactionNames, ReactionFormulas,'grRuleList',GeneNames);
%% Show N, metabolites, reactions, lower bounds, upper bounds

model=orig_model;
disp('S'); full(model.S)
disp('metabolites'); model.mets
disp('reactions'); model.rxns'
disp('lower bounds'); model.lb'
disp('upper bounds'); model.ub'
%% printing out reaction formulas

printRxnFormula(model,model.rxns);
%% set boundary condiations, e.g. measured fluxes

model = changeRxnBounds(model, {'R1'}, [1], 'b');
model = changeRxnBounds(model, {'R7'}, [0], 'l');
%% b) FBA (Flux Balance Analysis) simulating maximal growth

model=changeObjective(model,'mu');
FBAsolution=optimizeCbModel(model);
printFluxVector(model, FBAsolution.x, false, false)
[Involved_mets, Dead_ends] = draw_by_rxn(model, model.rxns, 'true', 'struc', {''}, {''}, FBAsolution.x)
%% c) Flux variability

[minFlux1,maxFlux1] = fluxVariability(model,100);
rxnNames = model.rxns; %{'PGI','PFK','FBP','FBA','TPI','GAPD','PGK','PGM','ENO','PYK','PPS','G6PDH2r','PGL','GND','RPI','RPE','TKT1','TKT2','TALA'};
rxnID = findRxnIDs(model,rxnNames);
printLabeledData(model.rxns(rxnID),[minFlux1(rxnID) maxFlux1(rxnID) maxFlux1(rxnID)-minFlux1(rxnID)],false)
%% Flux variability with R3=0

model = changeRxnBounds(model, {'R3','R6'}, [0; -1000], 'l');
model = changeRxnBounds(model, {'R3','R6'}, [0; 1000], 'u');
[minFlux,maxFlux] = fluxVariability(model,100);
rxnNames = model.rxns; %{'PGI','PFK','FBP','FBA','TPI','GAPD','PGK','PGM','ENO','PYK','PPS','G6PDH2r','PGL','GND','RPI','RPE','TKT1','TKT2','TALA'};
rxnID = findRxnIDs(model,rxnNames);
printLabeledData(model.rxns(rxnID),[minFlux(rxnID) maxFlux(rxnID) maxFlux(rxnID)-minFlux(rxnID)],false)
%% Flux variability with R3,R6=0

model = changeRxnBounds(model, {'R3','R6'}, [0; 0], 'l');
model = changeRxnBounds(model, {'R3','R6'}, [0; 0], 'u');
[minFlux,maxFlux] = fluxVariability(model,100);
rxnNames = model.rxns; %{'PGI','PFK','FBP','FBA','TPI','GAPD','PGK','PGM','ENO','PYK','PPS','G6PDH2r','PGL','GND','RPI','RPE','TKT1','TKT2','TALA'};
rxnID = findRxnIDs(model,rxnNames);
printLabeledData(model.rxns(rxnID),[minFlux(rxnID) maxFlux(rxnID) maxFlux(rxnID)-minFlux(rxnID)],false)
%% d) Analyzing flux correlations

model=orig_model;
model=changeRxnBounds(model, {'R1','R7'}, [1; 0], 'l');
model=changeRxnBounds(model, {'R1'}, [1], 'u');
model=changeObjective(model,'mu',[1]);
sol = optimizeCbModel(model);
growthRate = sol.f;
model = changeRxnBounds(model,'mu',0.9*growthRate,'l');

options.nStepsPerPoint = 8 * size(model.S, 2);
options.nPointsReturned = 1000;
  
[modelSampling1,samples1] = sampleCbModel(model,'smallnetwork_flux','ACHR',options);
rxnNames = {'mu','RBC1','RBC2','R2','R3','R4','R5','R6','R7'};
figure
sampleScatterMatrix(rxnNames,modelSampling1,samples1);
%% Visualization of the sampling results (correlation between reactions)

figure
for counter=1:size(samples1,1)
    for counter2=counter:size(samples1,1)
        if counter~=counter2
        subplot(size(samples1,1),size(samples1,1),(counter-1)*size(samples1,1)+counter2)
        plot(samples1(counter,:),samples1(counter2,:),'*')
        xlabel(modelSampling1.rxns(counter)),ylabel(modelSampling1.rxns(counter2))
        end
    end
end
%% Visualization of the distribution of the sampling results

plotSampleHist(modelSampling1.rxns, {samples1}, modelSampling1,[],[3,3]);
%% Identify sets of correlated reactions from the sampling data.

[sets,setNumber,setSize]=identifyCorrelSets(modelSampling1,samples1,1-1e-8);
for counter=1:size(sets)
    temp=cell2mat(sets(counter));
    disp(['Set ' num2str(counter) ':'])
    disp([num2cell(temp.set) temp.names])
end
%% What happens if you fix R3=0?

model = changeRxnBounds(model, {'R3'}, [0], 'b'); % fix R3=0
 
[modelSampling,samples] = sampleCbModel(model,'smallnetwork_flux','ACHR',options);
rxnNames = {'mu','RBC1','RBC2','R2','R4','R5','R6','R7'};
figure
sampleScatterMatrix(rxnNames,modelSampling,samples);
%% Visualization of the sampling results (correlation between reactions)

figure
for counter=1:size(samples,1)
    for counter2=counter:size(samples,1)
        if counter~=counter2
        subplot(size(samples,1),size(samples,1),(counter-1)*size(samples,1)+counter2)
        plot(samples(counter,:),samples(counter2,:),'*')
        xlabel(modelSampling.rxns(counter)),ylabel(modelSampling.rxns(counter2))
        end
    end
end
%% Visualization of the distribution of the sampling results

plotSampleHist(modelSampling.rxns, {samples}, modelSampling,[],[3,3]);
%% Identify sets of correlated reactions from the sampling data.

[sets,setNumber,setSize]=identifyCorrelSets(modelSampling,samples,1-1e-8);
for counter=1:size(sets)
    temp=cell2mat(sets(counter));
    disp(['Set ' num2str(counter) ':'])
    disp([num2cell(temp.set) temp.names])
end
%% What happens if you fi?x R3 = R6 = 0?

model = changeRxnBounds(model, {'R6'}, [0], 'b'); % fix R3=R6=0
  
[modelSampling,samples] = sampleCbModel(model,'smallnetwork_flux','ACHR',options);
rxnNames = {'mu','RBC1','RBC2','R2','R4','R5','R7'};
figure
sampleScatterMatrix(rxnNames,modelSampling,samples);
%% Visualization of the sampling results (correlation between reactions)

figure
for counter=1:size(samples,1)
    for counter2=counter:size(samples,1)
        if counter~=counter2
        subplot(size(samples,1),size(samples,1),(counter-1)*size(samples,1)+counter2)
        plot(samples(counter,:),samples(counter2,:),'*')
        xlabel(modelSampling.rxns(counter)),ylabel(modelSampling.rxns(counter2))
        end
    end
end
%% Visualization of the distribution of the sampling results

plotSampleHist(modelSampling.rxns, {samples}, modelSampling,[],[4,3]);
%% Identify sets of correlated reactions from the sampling data.

[sets,setNumber,setSize]=identifyCorrelSets(modelSampling,samples,1-1e-8);
for counter=1:size(sets)
    temp=cell2mat(sets(counter));
    disp(['Set ' num2str(counter) ':'])
    disp([num2cell(temp.set) temp.names])
end
%% e) Represent on the same plot the results from the FBA, FVA and sampling using R1 = 1 and R7 ? 0
% for the followings reactions: GR, R1, R2, R4, R5 and R6.

[minFlux1,maxFlux1] = fluxVariability(model,90);

rxns={'mu','R1','R2','R4','R5'};
figure;
for i = 1 : numel(rxns)
 subplot(2,3,i)
 pos=find(ismember(modelSampling1.rxns,rxns(i)));
%  h=histogram(samples1(pos,:),40);
%  h.FaceColor = 'k';
%  h.EdgeColor = 'k';
%  h.FaceAlpha=1;
[counts,centers]=hist(samples1(pos,:),15);
 plot(centers,counts,'k','LineWidth',2)
 hold on
 pos=find(ismember(model.rxns,rxns(i)));
 fig_axis=get(gca);
 max_y=fig_axis.YLim(2);
 plot([minFlux1(pos) minFlux1(pos)], [0 max_y],'r','LineWidth',2);
 plot([maxFlux1(pos) maxFlux1(pos)], [0 max_y],'r','LineWidth',2);
 plot([FBAsolution.x(pos) FBAsolution.x(pos)], [0 max_y],'y*','LineWidth',2);
 title(model.rxns{pos});
 xlabel('Flux')
 ylabel('# samples')
 hold off
end
legend('Sampling','MinFVA','MaxFVA','FBA')
%% f) Perform a single gene deletion study with the objective function mu. Discuss the result.

model=orig_model;
model=changeRxnBounds(model, {'R1','R7'}, [1; 0], 'l');
model=changeRxnBounds(model, {'R1'}, [1], 'u');
model = changeObjective(model,{'mu'},1);
% Single gene deletion study
[grRatio,grRateKO,grRateWT] = singleGeneDeletion(model,'FBA');
grRatio
%% g) Perform a double gene deletion study. What is the only non-lethal double deletion?
% Double gene deletion study

[grRatio,grRateKO,grRateWT] = doubleGeneDeletion(model,'FBA');
disp(['_' model.genes'; model.genes num2cell(grRatio)])