clear 
rxns_to_link='Lutein';
%In this script we want to enable the model to produce Luteine and S_epsilon_Carotene
% two demand reaction were added : epsilon carotene harvesting Luteinharvesting
epsilon=1e-4;
%set the path
path='C:\Users\maria.pacheco\Documents\GitHub\basicAnalysis';
%path='C:\Users\maria.pacheco\OneDrive - University of Luxembourg\Documents\GitHub\basicAnalysis';
%fix the names to facilitate readibility and fix compartemets that should
%be one letter within square brackets
model=readCbModel(strcat(path,'\Model\AraGEMconstraintFree.xml'));
model.mets= strrep(model.mets,'S__40_','S-');
model.mets= strrep(model.mets,'_44_','-');

model.mets= strrep(model.mets,'_41__45_','-');
model.mets= strrep(model.mets,'_45_','-');
model.mets= strrep(model.mets,'_c[C_c]','[c]');
model.mets= strrep(model.mets,'_p[C_p]','[p]');
model.mets= strrep(model.mets,'_m[C_m]','[m]');
model.mets= strrep(model.mets,'_x[C_x]','[x]');
model.mets= strrep(model.mets,'_v[C_v]','[v]');
model.mets= strrep(model.mets,'_ext[C_ext]','[e]');
model.mets= strrep(model.mets,'_acc[C_acc]','[a]');

model.mets= strrep(model.mets,'_biomass[C_biomass]','[b]');


%load(strcat(path,'\Model\Arabidopsis.mat'), 'model');
[TableNumerics, Required4biomass, Results, minimalMedia, model, consistentModel1]=InitialQualityControl(model, 'BIO_L');
disp(TableNumerics)

changeCobraSolver('ibm_cplex');
% check if the model can produce Lutein
rxns=findRxnsFromMets(model,model.mets(contains(model.metNames, rxns_to_link)));
disp(rxns)
printRxnFormula(model, rxns);
model=addReaction(model, 'Luteinharvesting',  'S_Lutein[c] ->');
model=changeObjective(model,'Luteinharvesting');
sol_Lut=optimizeCbModel(model,'max','zero');
disp(sol_Lut.f) %0

%the model can produce biomass but no lutein

% add demands to deads
if sol_Lut.f < epsilon
%add sink reaction to allow a flux and see where the networks has gaps
[ model_corrected,deadEnd]=findDeadEndsFastbox(model);

model=model_corrected;
clear model_corrected
model=changeObjective(model,'Luteinharvesting');
sol_Lut=optimizeCbModel(model,'max','zero');
disp(sol_Lut.f) %500

model=changeObjective(model,'BIO_L');
sol_Bio=optimizeCbModel(model,'max','zero');
disp(sol_Bio.f)%  433
consistentModel1=changeObjective(consistentModel1,'BIO_L');
sol_BioCons=optimizeCbModel(consistentModel1,'max','zero');
disp(sol_BioCons.f) % 134

rxnsRequired4Obj=model.rxns(sol_Lut.x~=0);
Table1=table(rxnsRequired4Obj, printRxnFormula(model,rxnsRequired4Obj),sol_Lut.x(ismember(model.rxns,rxnsRequired4Obj)));

% Table 1 shows that we form Lutein from Isopentenyl di p which correspond
% to the literature 'S_Dehydrodolichol_32_diphosphate_c[C_c]  <=> ' that
% enters from a sink is being converted into Isopentenyl

% so we close this input

model.lb(ismember(model.rxns,'SinkToRemoveDeadEnd_I__S_Dehydrodolichol_32_diphosphate[c]'))=0;
model=changeObjective(model,'Luteinharvesting');
sol_Lut=optimizeCbModel(model,'max','zero');
disp(sol_Lut.f)% 500

model=changeObjective(model,'BIO_L');
sol_Bio=optimizeCbModel(model,'max','zero');
disp(sol_Bio.f) %428


%By closing 

rxnsRequired4Obj=model.rxns(sol_Lut.x~=0);
Table2=table(rxnsRequired4Obj, printRxnFormula(model,rxnsRequired4Obj),sol_Lut.x(ismember(model.rxns,rxnsRequired4Obj)));
model.lb(ismember(model.rxns,'SinkToRemoveDeadEnd_I__S_2-C-Methyl-D-erythritol_32_2-4-cyclodiphosphate[c]'))=0;
model=changeObjective(model,'Luteinharvesting');
sol_Lut=optimizeCbModel(model,'max','zero');
disp(sol_Lut.f)%250
rxns=findRxnsFromMets(model,'S-R-5-Diphosphomevalonate[c]');
printRxnFormula(model, rxns)
% Connection to Mevalonate is missing
%add reaction to connect to diphosphomevalonate

model=addReaction(model, 'diphosphomevalonate decarboxylase', 'S_ATP[c]  + S-R-5-Diphosphomevalonate[c] -> S_ADP[c] + S_Isopentenyl_32_diphosphate[c] + S_Pyrophosphate[c] +  S_CO2[m]');
model.lb(ismember(model.rxns,'SinkToRemoveDeadEnd_I__S-R-5-Diphosphomevalonate[c]'))=0;
model=changeObjective(model,'Luteinharvesting');
sol_Lut=optimizeCbModel(model,'max','zero');
disp(sol_Lut.f) %250

model=changeObjective(model,'BIO_L');
sol_Bio=optimizeCbModel(model,'max','zero');
disp(sol_Bio.f) %428

rxnsRequired4Obj=model.rxns(sol_Lut.x~=0);
Table3=table(rxnsRequired4Obj, printRxnFormula(model,rxnsRequired4Obj),sol_Lut.x(ismember(model.rxns,rxnsRequired4Obj)));
Table3(contains(Table3.Var2,'Meva'),:)
Table3(contains(Table3.Var2,'S-S-3-Hydroxy-3-methylglutaryl-CoA'),:)

% the lutein is connected to the mevalonate pathway and the pathways starts
% at AcetylCoA

%Pyruvate is converted into acetylcoa in the mitocondria however no reactions 
%transport accetycoa in the cytoplasma 

rxns=findRxnsFromMets(model,'S_Pyruvate[m]');
printRxnFormula(model, rxns)
rxns2=findRxnsFromMets(model,'S_Acetyl-CoA[m]');
printRxnFormula(model, rxns2)


% lets add it

model=addReaction(model,'S_Acetyl_45_CoA_trans', 'S_Acetyl-CoA[m] <=> S_Acetyl-CoA[c]');
model=changeObjective(model,'Luteinharvesting');
sol_Lut=optimizeCbModel(model,'max','zero');
disp(sol_Lut.f)%250


rxnsRequired4Obj=model.rxns(sol_Lut.x~=0);
Table4=table(rxnsRequired4Obj, printRxnFormula(model,rxnsRequired4Obj),sol_Lut.x(ismember(model.rxns,rxnsRequired4Obj)));
%We should close non necessary sinks
model=changeObjective(model,'BIO_L');
sol_Bio=optimizeCbModel(model,'max','zero');
disp(sol_Bio.f) 

model.lb(contains(model.rxns,'SinkToRemove'))=0;
model.lb(contains(model.rxns,'dead_end_II'))=0;

model=changeObjective(model,'Luteinharvesting');
sol_Lut=optimizeCbModel(model,'max','zero');
disp(sol_Lut.f)

rxnsRequired4Obj=model.rxns(sol_Lut.x~=0);
Table5=table(rxnsRequired4Obj, printRxnFormula(model,rxnsRequired4Obj),sol_Lut.x(ismember(model.rxns,rxnsRequired4Obj)));

% the flux is comminng fom glutamine and glutamate
% close glutamine and glutamine 
%Lets close to see where flux can come from
tmp=model.lb(contains(model.rxns,'Ex6'));
tmp2=model.lb(contains(model.rxns,'Ex7'));
tmp3=model.lb(contains(model.rxns,'Ex9'));
tmp4=model.lb(contains(model.rxns,'Ex10'));


model.lb(contains(model.rxns,'Ex6'))=0;
model.lb(contains(model.rxns,'Ex7'))=0;
model.lb(contains(model.rxns,'Ex9'))=0;
model.lb(contains(model.rxns,'Ex10'))=0;

sol_Lut=optimizeCbModel(model,'max','zero');

sol_Lut.f
rxnsRequired4Obj=model.rxns(sol_Lut.x~=0);
Table6=table(rxnsRequired4Obj, printRxnFormula(model,rxnsRequired4Obj),sol_Lut.x(ismember(model.rxns,rxnsRequired4Obj)));
% all the carbon inputs are closed but there is still fluxes showing that
% there is sth wrong with the model
%The Flux is coming F6P but then is spinning in a loop
rxns=findRxnsFromMets(model,'S_beta-D-Fructose_32_6-phosphate[c]');
printRxnFormula(model,rxns)

rxns=findRxnsFromMets(model,'S_alpha-D-Glucose_32_6-phosphate[c]');
printRxnFormula(model,rxns);
rxns=findRxnsFromMets(model,'S_alpha-D-Glucose[c]');
printRxnFormula(model,rxns);
 % put Glutamate and Glutamine back in
model.lb(contains(model.rxns,'Ex6'))=tmp;
model.lb(contains(model.rxns,'Ex7'))=tmp2;
model.lb(contains(model.rxns,'Ex9'))=tmp3;
model.lb(contains(model.rxns,'Ex10'))=tmp4;

model=changeObjective(model,'BIO_L');
sol_Bio=optimizeCbModel(model,'max','zero');
disp(sol_Bio.f )
% Lets allow the uptake of alpha glucose
model.lb(contains(model.rxns,'Ex13'))=-1000;
model=changeObjective(model,'Luteinharvesting');

sol_Lut=optimizeCbModel(model,'max','zero');

disp(sol_Lut.f)
model=changeObjective(model,'BIO_L');
sol_Bio=optimizeCbModel(model,'max','zero');
sol_Bio.f 

rxnsRequired4Obj=model.rxns(sol_Lut.x~=0);
Table7=table(rxnsRequired4Obj, printRxnFormula(model,rxnsRequired4Obj),sol_Lut.x(ismember(model.rxns,rxnsRequired4Obj)));

%Lets remove some sinks and connectings
rxns=model.rxns(contains(model.rxns,'SinkToRemoveDeadEnd'));
rxns=setdiff(rxns, Table7.rxnsRequired4Obj);
model=removeRxns(model,rxns);

rxns=findRxnsFromMets(model,'S_CO2[m]');
printRxnFormula(model,rxns)
model.mets(contains(model.mets, 'CO2'))
% lets allow CO2 to be exported and imported

model=addReaction(model,'S_CO2_trans_m_c', 'S_CO2[m] <=> S_CO2[c]');
model=addReaction(model,'S_CO2_trans_c_e', 'S_CO2[c] <=> S_CO2[e]');

model=addReaction(model,'Ex_CO2', 'S_CO2[e] <=> ');

model=removeRxns(model,'SinkToRemoveDeadEnd_I__S_CO2[m]');
model=removeRxns(model,'SinkToRemoveDeadEnd_I__S_CO2[p]');

rxns=findRxnsFromMets(model,'S_Hydrogen[c]');
printRxnFormula(model,rxns)
model.mets(contains(model.mets, 'CO2'))
% lets allow CO2 to be exported and imported
model=addReaction(model,'S_H_trans_m_e', 'S_Hydrogen[c] <=> S_Hydrogen[e]');
model=addReaction(model,'Ex_H', 'S_Hydrogen[e] <=> ');
model=removeRxns(model,'dead_end_II_S_Hydrogen[c]');
model=changeObjective(model,'Luteinharvesting');
sol_Lut=optimizeCbModel(model,'max','zero');
disp(sol_Lut.f)
model=changeObjective(model,'BIO_L');
sol_Bio=optimizeCbModel(model,'max','zero');
disp(sol_Bio.f)
rxnsRequired4Obj=model.rxns(sol_Lut.x~=0);
Table8=table(rxnsRequired4Obj, printRxnFormula(model,rxnsRequired4Obj),sol_Lut.x(ismember(model.rxns,rxnsRequired4Obj)));
end
[TableNumerics2, Required4biomass2, Results2, minimalMedia2, model2, consistentModel2]=InitialQualityControl(model, 'BIO_L');
save manualGapfilling


