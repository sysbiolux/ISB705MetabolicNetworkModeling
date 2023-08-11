%driver script
clear 
addpath(genpath('C:\Program Files\IBM\ILOG\\CPLEX_Studio1210\cplex\matlab\x64_win64'))
changeCobraSolver('ibm_cplex')
%model=readCbModel('C:\Users\maria.pacheco\OneDrive - University of Luxembourg\Documents\GitHub\ISB705MetabolicNetworkModelling\Model\Blautia_hydrogenotrophica_DSM_10507.xml');

model=readCbModel('C:\Users\maria.pacheco\Documents\GitHub\basicAnalysis\Model\Blautia_hydrogenotrophica_DSM_10507.xml');
[TableNumerics, Required4biomass, Results, minimalMedia, model, consistentModel]=InitialQualityControl(model, 'EX_biomass(e)');


load('C:\Users\maria.pacheco\Documents\GitHub\basicAnalysis\Model\brain_model.mat')
[TableNumerics2, Required4biomass2, Results2, minimalMedia2, model2]=InitialQualityControl(brain_model, 'biomass_maintenance');
