clear 

%In this script we want to enhable the model to produce Luteine and S_epsilon_45_Carotene
% two demand reaction were added : epsilon carotene harvesting Luteinharvesting
epsilon=1e-4;
%set the path
path='C:\Users\maria.pacheco\OneDrive - University of Luxembourg\Documents\GitHub\basicAnalysis';
model=readCbModel(strcat(path,'\Model\AraGEMconstraintFree.xml'));
%load(strcat(path,'\Model\Arabidopsis.mat'), 'model');
changeCobraSolver('ibm_cplex');
% check if the model can produce Lutein

[TableNumerics, Required4biomass, Results, minimalMedia, model]=InitialQualityControl(model, 'BIO_L');
rxns=findRxnsFromMets(model,model.mets(contains(model.metNames, 'Lutein')));
disp(rxns)
model=addReaction(model, 'Starch_stock',  'S_Starch_p[C_p] <=>');

model=addReaction(model, 'Luteinharvesting',  'S_Lutein_c[C_c] ->');
model=changeObjective(model,'Luteinharvesting');
sol_ori=optimizeCbModel(model,'max','zero');
sol_ori.f
% add demands to deadend
if sol_ori.f < epsilon
[ model_corrected,deadEnd]=findDeadEndsFastbox(model);
% model_corrected.rxnNames(end+1:numel(model_corrected.rxns))=model_corrected.rxns(numel(model_corrected.rxnNames)+1:numel(model_corrected.rxns));
% model_corrected.c(end+1:numel(model_corrected.rxns))=0;
model=model_corrected;
sol_ori=optimizeCbModel(model,'max','zero');
sol_ori.f
end
rxnsRequired4Obj=model.rxns(sol_ori.x~=0);
T=table(rxnsRequired4Obj, printRxnFormula(model,rxnsRequired4Obj),sol_ori.x(ismember(model.rxns,rxnsRequired4Obj)));
% table 1 shows that we form Lutein from Isopentenyl di p which correspond
% to the literature 'S_Dehydrodolichol_32_diphosphate_c[C_c]  <=> ' that
% enters from a sink is being converted into Isopentenyl
% so we close this input
model.lb(ismember(model.rxns,'SinkToRemoveDeadEnd_I_S_G13032_c[C_c]_S_G13032_c[C_c]'))=0;
sol_ori=optimizeCbModel(model,'max','zero');
sol_ori.f
rxnsRequired4Obj=model.rxns(sol_ori.x~=0);
T2=table(rxnsRequired4Obj, printRxnFormula(model,rxnsRequired4Obj),sol_ori.x(ismember(model.rxns,rxnsRequired4Obj)));
% the lutein is connected to the mevalonate pathway and the pathways starts
% at AcetylCoA
%Now we need to connect it to the glycolyse because acetylCoa comes from
%'S_Benzoyl_32_acetyl_45_CoA_c[C_c]  <=> '
%Pyruvate is converted into acetylcoa in the mitocondria however no reactions 
%transport accetycoa in the cytoplasma 
rxns=findRxnsFromMets(model,'S_Pyruvate_m[C_m]');
printRxnFormula(model, rxns)
% lets add it
model=addReaction(model,'S_Acetyl_45_CoA_trans', 'S_Acetyl_45_CoA_c[C_c] <=> S_Acetyl_45_CoA_m[C_m]');

%lets force the flux to come from pyruvate
sol_ori=optimizeCbModel(model,'max','zero');

sol_ori.f
rxnsRequired4Obj=model.rxns(sol_ori.x~=0);
T3=table(rxnsRequired4Obj, printRxnFormula(model,rxnsRequired4Obj),sol_ori.x(ismember(model.rxns,rxnsRequired4Obj)));
%By adding sink reactions we allowed fluxes to enter the systems we need to
%close the input that are not needed
model_org=model;
val=T3.Var3(ismember(T3.rxnsRequired4Obj,'R01978_c'));
model.lb(ismember(model.rxns,'R00209_m'))=abs(val);
model.lb(ismember(model.rxns,'R00200_c'))=abs(val);

model.lb(end)=- abs(val);
model.ub(end)=-abs(val);

%Lets force all the flux to go through the transporter we added

sol_ori=optimizeCbModel(model,'max','zero');

sol_ori.f
rxnsRequired4Obj=model.rxns(sol_ori.x~=0);
T4=table(rxnsRequired4Obj, printRxnFormula(model,rxnsRequired4Obj),sol_ori.x(ismember(model.rxns,rxnsRequired4Obj)));
% now the flux is coming from Glucose but we should close non necessary sinks
model=model_org;
model.lb(contains(model.rxns,'SinkToRemove'))=0;
model.lb(contains(model.rxns,'dead_end_II'))=0;

sol_ori=optimizeCbModel(model,'max','zero');

sol_ori.f;
rxnsRequired4Obj=model.rxns(sol_ori.x~=0);
T5=table(rxnsRequired4Obj, printRxnFormula(model,rxnsRequired4Obj),sol_ori.x(ismember(model.rxns,rxnsRequired4Obj)));


